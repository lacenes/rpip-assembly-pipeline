#!/usr/bin/env python3

"""
Pipeline Python para montagem gene-a-gene
-------------------------------------------------------
Passos principais:
1. Extrai cada gene de um multifasta (samtools faidx)
2. Indexa o gene (BWA) e alinha as leituras emparelhadas
3. Gera BAM ordenado + idxstats + depth
4. Chama variantes (bcftools) e cria FASTA de consenso (opcional máscara < min_cov)
5. Recruta apenas as leituras que alinham (FLAG 0x2) → FASTQ pareado
6. Monta com metaSPAdes
7. Consolida:
   • all_spades_scaffolds.fasta   (headers: <gene>_scaffoldN)
   • all_consensus_genes.fasta    (headers: <gene>_consensus)
8. Reverse-complement para o multifasta gerado pelo spades
"""
import subprocess
import shutil
import sys
import logging
import argparse
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

try:
    import argcomplete
except ImportError:
    argcomplete = None


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )


def check_tools(tools):
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        logging.error("Ferramentas faltando no PATH: %s", ", ".join(missing))
        sys.exit(1)


def run(cmd, **kwargs):
    logging.info("Running: %s", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True, **kwargs)
    except subprocess.CalledProcessError as e:
        logging.error("Falha ao executar: %s", e.cmd)
        sys.exit(1)


def reverse_complement(seq: str) -> str:
    complement = {'A':'T','T':'A','C':'G','G':'C',
                  'a':'t','t':'a','c':'g','g':'c',
                  'N':'N','n':'n'}
    return ''.join(complement.get(b, 'N') for b in reversed(seq))


def process_reverse_fasta(input_fasta: Path, output_fasta: Path) -> None:
    """
    Gera reverse-complement de todos os scaffolds.
    """
    with open(output_fasta, 'w') as out_f:
        for record in SeqIO.parse(str(input_fasta), 'fasta'):
            rid = record.id
            if '_scaffold' in rid:
                rc_seq = reverse_complement(str(record.seq))
                rc_record = record.__class__(
                    id=rid + '_rc',
                    seq=Seq(rc_seq),
                    description=''
                )
                SeqIO.write(rc_record, out_f, 'fasta')


def main():
    setup_logging()
    check_tools(['samtools','bwa','bcftools','spades.py'])

    parser = argparse.ArgumentParser("Pipeline gene-a-gene (metagenoma)")
    parser.add_argument("-g", "--genes-fasta", required=True, help="Multifasta com genes‐alvo")
    parser.add_argument("-p", "--prefix-base", required=True, help="Prefixo base para saídas")
    parser.add_argument("-1", "--r1", dest="r1", required=True, help="FASTQ R1")
    parser.add_argument("-2", "--r2", dest="r2", required=True, help="FASTQ R2")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Threads")
    parser.add_argument("--min-cov", type=int, default=10, help="Profundidade mínima para consenso")
    if argcomplete:
        argcomplete.autocomplete(parser)
    args = parser.parse_args()

    genes_fa = Path(args.genes_fasta)
    prefix_b = args.prefix_base
    r1 = args.r1
    r2 = args.r2
    threads = str(args.threads)
    min_cov = args.min_cov

    outdir = Path("run_by_gene")
    outdir.mkdir(exist_ok=True)

    # 0) Index multifasta
    fai = genes_fa.with_suffix(genes_fa.suffix + ".fai")
    if not fai.exists():
        run(["samtools","faidx",str(genes_fa)])

    genes = [ln.split("\t")[0] for ln in open(fai)]

    # Loop por gene
    for gene in genes:
        logging.info("=== Gene %s ===", gene)
        gdir = outdir / gene
        gdir.mkdir(exist_ok=True)
        gfa = gdir / f"{gene}.fasta"
        prefix = gdir / f"{prefix_b}_{gene}"

        # 1) Extrair gene
        run(["samtools","faidx",str(genes_fa),gene,"-o",str(gfa)])

        # 2) Index BWA
        if not gfa.with_suffix(gfa.suffix + ".bwt").exists():
            run(["bwa","index",str(gfa)])

        # 3) Alinhamento → BAM ordenado
        bam = prefix.with_suffix('.sorted.bam')
        p1 = subprocess.Popen([
            "bwa","mem","-t",threads,str(gfa),r1,r2
        ], stdout=subprocess.PIPE)
        p2 = subprocess.Popen([
            "samtools","view","-bS","-@",threads,"-"
        ], stdin=p1.stdout, stdout=subprocess.PIPE)
        run(["samtools","sort","-@",threads,"-o",str(bam)], stdin=p2.stdout)
        run(["samtools","index",str(bam)])

        # 4) Estatísticas + depth
        run(["samtools","idxstats",str(bam)], stdout=open(prefix.with_suffix('.idxstats.tsv'),'w'))
        depth_file = prefix.with_suffix('.depth.txt')
        run(["samtools","depth","-a",str(bam)], stdout=open(depth_file,'w'))

        # 5) Variantes + consenso
        vcf = prefix.with_suffix('.variants.vcf.gz')
        p3 = subprocess.Popen([
            "bcftools","mpileup","-f",str(gfa),"-Ou",str(bam)
        ], stdout=subprocess.PIPE)
        run(["bcftools","call","--ploidy","1","-m","-v","-Oz","-o",str(vcf)], stdin=p3.stdout)
        run(["bcftools","index",str(vcf)])

        # 6) Máscara de baixa cobertura
        bed = prefix.with_suffix('.lowcov.bed')
        with open(depth_file) as dfh, open(bed,'w') as bedfh:
            for ln in dfh:
                ctg,pos,cov = ln.split()[:3]
                if int(cov) < min_cov:
                    bedfh.write(f"{ctg}\t{int(pos)-1}\t{pos}\n")
        mask_opt = ["-m",str(bed)] if min_cov else []

        # 7) Gera consenso
        consensus = gdir / f"{prefix_b}_{gene}_consensus.fasta"
        run(["bcftools","consensus","-f",str(gfa)] + mask_opt + [str(vcf)], stdout=open(consensus,'w'))

        # 8) Extrai pares mapeados e SPAdes
        paired = prefix.with_suffix('.paired.bam')
        namesort = prefix.with_suffix('.namesorted.bam')
        r1_out = gdir / f"{prefix.name}_R1.fastq"
        r2_out = gdir / f"{prefix.name}_R2.fastq"
        run(["samtools","view","-b","-f","0x2",str(bam),"-o",str(paired)])
        run(["samtools","sort","-n","-@",threads,str(paired),"-o",str(namesort)])
        run(["samtools","fastq","-1",str(r1_out),"-2",str(r2_out),"-0","/dev/null","-s","/dev/null","-n",str(namesort)])

        spdir = gdir / f"{prefix_b}_{gene}_spades"
        try:
            logging.info("Running SPAdes for gene %s", gene)
            subprocess.run([
                "spades.py","--meta",
                "-k","21,33,37,39",
                "-1",str(r1_out),
                "-2",str(r2_out),
                "-o",str(spdir),
                "--threads",threads,
                "-m","64"
            ], check=True)
        except subprocess.CalledProcessError as e:
            logging.error("SPAdes failed for gene %s: %s", gene, e)
            continue

    # 9) Consolida scaffolds & consensos
    all_scaf = outdir / "all_spades_scaffolds.fasta"
    all_cons = outdir / "all_consensus_genes.fasta"
    with open(all_scaf,'w') as scaf_out, open(all_cons,'w') as cons_out:
        for gene in genes:
            sfile = outdir / gene / f"{prefix_b}_{gene}_spades" / "scaffolds.fasta"
            if sfile.exists():
                n = 0
                for ln in open(sfile):
                    if ln.startswith('>'):
                        n += 1
                        scaf_out.write(f">{gene}_scaffold{n}\n")
                    else:
                        scaf_out.write(ln)
            cfile = outdir / gene / f"{prefix_b}_{gene}_consensus.fasta"
            if cfile.exists():
                for ln in open(cfile):
                    if ln.startswith('>'):
                        cons_out.write(f">{gene}_consensus\n")
                    else:
                        cons_out.write(ln)

    # 10) Reverse-complement dos scaffolds marcados '_rv'
    all_scaf_rc = outdir / "all_spades_scaffolds_rc.fasta"
    process_reverse_fasta(all_scaf, all_scaf_rc)

    logging.info("Pipeline concluído. Arquivos principais:")
    logging.info(" • %s", all_scaf)
    logging.info(" • %s", all_scaf_rc)
    logging.info(" • %s", all_cons)


if __name__ == "__main__":
    main()
