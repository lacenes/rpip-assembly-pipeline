# Pipeline de Montagem Gene-a-Gene utilizado nos dados brutos do RPIP (Illumina)

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Status](https://img.shields.io/badge/status-stable-brightgreen)]()

## 📌 Visão Geral
Pipeline para montagem direcionada (gene-a-gene) a partir de leituras pareadas.  
Este pipeline automatiza a montagem de genes específicos a partir de dados obtidos por painéis de enriquecimento por captura ou por metagenômica, utilizando um multifasta de genes-alvo como referência.

### **Etapas principais:**
- Extração de genes-alvo a partir do multifasta (Samtools)  
- Indexação e alinhamento das leituras para cada gene (BWA)  
- Geração de BAM ordenado, estatísticas de cobertura e profundidade (Samtools)  
- Chamada de variantes e geração de sequência consenso com máscara para regiões de baixa cobertura (BCFtools) > montagem por referência
- Recrutamento de leituras pareadas mapeadas para montagem direcionada (metaSPAdes)   > montagem de novo
- Consolidação de scaffolds e consensos, incluindo geração de reverse-complement dos scaffolds do SPAdes  

---

### ⚙️ Requisitos
- Python 3.8+
- [Samtools](http://www.htslib.org/)  
- [BWA](http://bio-bwa.sourceforge.net/)  
- [BCFtools](http://samtools.github.io/bcftools/)  
- [SPAdes](https://ablab.github.io/spades/installation.html)
- [Biopython](https://biopython.org/)
- [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main)
- [Seqkit](https://bioinf.shenwei.me/seqkit/)

Certifique-se de que todas as ferramentas estejam disponíveis no `$PATH`.

---

### 📥 Instalação

### Criação do ambiente Conda (recomendada)
```bash
conda env create -f assembly_env_minimal.yml
conda activate assembly
```

### Instalação manual de dependências Python (opcional)
```bash
pip install biopython argcomplete
```

### Clonagem do repositório
```bash
git clone https://github.com/lacenes/rpip-assembly-pipeline.git
cd rpip-assembly-pipeline
```

---

### ▶️ Uso
```bash
python assembly_by_gene_rc.py \
  -g genes_alvo.fasta \
  -p amostra1 \
  -1 amostra1_R1.fastq.gz \
  -2 amostra1_R2.fastq.gz \
  -t 18 \
  --min-cov 10
```

### Parâmetros
```text
  -g, --genes-fasta : Multifasta com genes-alvo
  -p, --prefix-base : Prefixo para arquivos de saída
  -1, --r1          : Leituras pareadas R1
  -2, --r2          : Leituras pareadas R2
  -t, --threads     : Número de threads (padrão: 8)
  --min-cov         : Cobertura mínima para consenso (padrão: 10)
```

---

### 📂 Estrutura de Saída
```text
run_by_gene/
   ├── <gene1>/
   │    ├── <prefixo>_<gene1>_consensus.fasta
   │    ├── <prefixo>_<gene1>_spades/
   │    └── ...
   ├── all_spades_scaffolds.fasta        (Montagem de novo) 
   ├── all_spades_scaffolds_rc.fasta     (Montagem de novo reverso complementar) 
   └── all_consensus_genes.fasta         (Montagem por referência)
```

---

### 📌 Citação

- **Softwares**: Cite todos os softwares usados nesse pipiline (Bloco Requisitos)
- **Autor**: Thiago de Jesus Sousa  
- **Email**: thiagojsousa@gmail.com  
- **LinkedIn**: [linkedin.com/in/thiago-sousa-12106b79](https://www.linkedin.com/in/thiago-sousa-12106b79)

---

### 📜 Licença
Distribuído sob a licença MIT. Consulte o arquivo [LICENSE](LICENSE) para mais informações.
