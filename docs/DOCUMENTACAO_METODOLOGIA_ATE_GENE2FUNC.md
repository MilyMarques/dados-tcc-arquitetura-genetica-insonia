# Metodologia do TCC até FUMA SNP2GENE, MAGMA e GENE2FUNC

**Projeto:** Caracterização da Arquitetura da Insônia  
**Curso:** Biomedicina, Unifesp/EPM  
**Escopo deste documento:** etapas concluídas desde a seleção dos estudos até a interpretação funcional inicial no FUMA GENE2FUNC. KEGG, Reactome, g:Profiler, ShinyGO e STRING ainda não fazem parte deste documento.

## 1. Objetivo analítico

O objetivo desta etapa foi processar de forma reproduzível estatísticas-resumo de GWAS de insônia, executar anotação e mapeamento gene-variável no FUMA SNP2GENE, obter associações no nível de genes por MAGMA e definir uma lista priorizada de genes para interpretação funcional posterior.

O estudo não interpreta uma variante GWAS como causal por si só. A atribuição de genes é tratada como priorização de evidência, combinando significância gene-baseada, mapeamento funcional e consistência entre estudos.

## 2. Estudos e papel de cada conjunto de dados

### Estudos principais

1. Hammerschlag et al. (2017), insônia.
2. Lane et al. (2019), fenótipo extremo de sintomas de insônia versus controles.
3. Jansen et al. (2019), componente publicamente disponível do UK Biobank.
4. Watanabe et al. (2022), componente publicamente disponível do UK Biobank.

Os resultados públicos de Jansen e Watanabe não contêm a parcela restrita de 23andMe. Eles foram mantidos porque são dados publicamente acessíveis e relevantes, mas sua sobreposição com UK Biobank foi considerada nas comparações.

### Camada de sensibilidade OpenGWAS

Foram também processados dois GWAS de insônia derivados de UK Biobank/OpenGWAS: `ukb-a-13` (Neale Lab) e `ukb-b-3957` (PHESANT). Eles foram usados para avaliar a consistência técnica da lista principal, sem acrescentar unidades independentes de evidência, pois representam amostras UK Biobank sobrepostas.

O conjunto FinnGen `finn-b-F5-INSOMNIA` foi processado e validado. Não apresentou SNPs com p < 5e-8 no arquivo disponibilizado; portanto, o FUMA não identificou SNPs candidatos. Este resultado foi registrado e o limiar não foi reduzido artificialmente.

## 3. Processamento das summary statistics

Os arquivos foram processados por scripts Python configuráveis, com execução por blocos para comportar arquivos GWAS grandes. Foram criados arquivos compactados `.txt.gz` específicos para FUMA e relatórios JSON de auditoria/QC.

Em cada estudo, o processamento incluiu:

- leitura do arquivo bruto e inspeção dos nomes de colunas;
- padronização das colunas de cromossomo, posição, alelo de efeito, alelo alternativo, p-valor, tamanho amostral e identificador de variante, quando disponível;
- conversão de cromossomos para o formato aceito pelo pipeline final (autossomos 1–22);
- garantia de que posições fossem inteiros em coordenadas hg19;
- exclusão de linhas sem p-valor válido, fora do intervalo 0 < p <= 1 ou sem dados obrigatórios;
- inspeção de duplicatas e de valores extremos;
- preservação de campos necessários para submissão;
- compactação real com gzip, e não apenas renomeação do arquivo.

O processo produziu um input FUMA para cada estudo e um relatório de validação. A validação verificou, entre outros itens, cabeçalho, delimitador, p-valores, coordenadas, duplicatas e colunas selecionadas no formulário do FUMA.

## 4. FUMA SNP2GENE

Cada um dos quatro estudos principais foi submetido separadamente ao SNP2GENE, com configurações harmonizadas.

### Identificação de loci

- Limiar de SNP independente significativo: p < 5e-8;
- limiar de SNP candidato: p < 0,05;
- r2 para SNPs independentes significativos: 0,6;
- r2 para lead SNPs: 0,1;
- painel de LD: 1000 Genomes Phase 3, população EUR;
- inclusão de variantes do painel de referência em LD: sim;
- distância máxima para unir blocos de LD: 250 kb.

### Mapeamento SNP2GENE

Foram preservadas as três modalidades de evidência fornecidas pelo FUMA:

1. **Mapeamento posicional:** associa variantes a genes pela localização genômica.
2. **Mapeamento eQTL:** relaciona variantes a genes cuja expressão é associada à variante; nas integrações, foi destacada a evidência cerebral.
3. **Mapeamento por interação cromatínica:** relaciona regiões em LD a promotores por dados Hi-C, com atenção a conjuntos corticais.

O mapeamento foi utilizado como evidência de priorização, e não como prova de causalidade. Na integração, somente contagens maiores que zero foram consideradas evidência real; valores zero nos arquivos do FUMA indicam ausência daquele tipo de mapeamento.

### Anotações funcionais e região MHC

Foram usadas anotações funcionais disponibilizadas pelo FUMA. A região MHC foi excluída apenas das anotações/mapeamentos, mantendo-a na análise MAGMA. Esta decisão evita que a estrutura de LD complexa do MHC domine a interpretação baseada em SNPs sem modificar o universo de testes de genes do MAGMA.

### Versão de genes

Os jobs principais de Jansen e Watanabe utilizaram Ensembl v110; o FinnGen corrigido e o GENE2FUNC do Core16 utilizaram Ensembl `v110_userGen`. Para o núcleo posterior, foram retidos genes *protein-coding*, garantindo uma lista interpretável nas ferramentas de vias e redes proteicas.

## 5. MAGMA

O MAGMA foi executado no FUMA para cada estudo. A significância por gene foi avaliada pelo p-valor do MAGMA e pela correção de Bonferroni específica ao número de genes testados em cada job. O FDR de Benjamini-Hochberg foi também calculado para tabelas descritivas.

Os resultados MAGMA foram usados como a camada estatística principal, pois testam associação gene-baseada em vez de apenas transferir a posição de um SNP ao gene próximo.

## 6. Integração entre estudos e definição do núcleo

Os resultados de lead SNPs, loci, genes mapeados e MAGMA foram extraídos dos quatro jobs do FUMA e padronizados em tabelas comparáveis. Foram comparados loci e genes entre estudos.

Para evitar dupla contagem da mesma base de participantes, Jansen e Watanabe foram agrupados como uma única unidade de evidência denominada `ukb_public`. Assim, havia três unidades principais: Hammerschlag, Lane e UKB público.

Um gene entrou no núcleo se atendesse simultaneamente aos seguintes critérios:

1. ser protein-coding;
2. ter significância MAGMA por Bonferroni em pelo menos duas unidades analíticas;
3. apresentar ao menos uma modalidade positiva de mapeamento SNP2GENE (posicional, eQTL cerebral ou interação cromatínica cortical).

Foram priorizados como alta prioridade os genes com mais de uma modalidade de mapeamento; genes com suporte MAGMA equivalente, mas menor diversidade de mapeamento, foram classificados como prioridade intermediária.

### Lista resultante

**Alta prioridade (12):** MEIS1, ACBD4, HEXIM2, HEXIM1, FOXP2, SMIM32, SMAD5, NMT1, DCAKD, PAM, KLC2 e DNPH1.

**Prioridade intermediária (4):** TCF4, PLCD3, TRPC7 e LIN28B.

O núcleo não deve ser descrito como uma lista de genes causais. Ele representa genes prioritários com evidência estatística e funcional convergente.

## 7. Sensibilidade com OpenGWAS UKB

Os dois jobs UKB OpenGWAS foram analisados separadamente. Os 16 genes do núcleo apresentaram significância por FDR nos dois conjuntos. Treze apresentaram significância por Bonferroni em ambos; PAM, KLC2 e LIN28B foram Bonferroni-significativos em uma das duas rodagens e FDR-significativos nas duas.

Essa etapa reforça a consistência técnica no contexto UK Biobank, mas não aumenta o número de estudos independentes e não modifica a lista núcleo.

## 8. GENE2FUNC

O GENE2FUNC foi executado com:

- lista de interesse: os 16 genes núcleo em identificadores Ensembl;
- background personalizado: 18.717 genes que foram efetivamente testáveis pelo MAGMA nos quatro estudos principais;
- Ensembl: `v110_userGen`;
- expressão: GTEx v8 (53 tecidos e 30 categorias gerais) e BrainSpan (idade e desenvolvimento);
- ajuste de múltiplos testes: Benjamini-Hochberg, FDR < 0,05;
- sobreposição mínima por conjunto: dois genes.

O background personalizado foi usado porque representa o universo de genes que poderiam ter sido escolhidos pelo procedimento MAGMA. Ele é metodologicamente preferível a comparar os 16 genes contra todo o genoma.

A opção MHC não foi reaplicada nesta etapa: GENE2FUNC recebe listas de genes, não loci genômicos, e alterar o background nessa fase quebraria a compatibilidade com o universo MAGMA.

### Resultado de expressão por tecido

Não houve enriquecimento de expressão diferencial significativo após FDR em GTEx ou BrainSpan. Isso deve ser interpretado como ausência de evidência de concentração tecidual detectável com a lista de 16 genes, e não como ausência de relevância cerebral ou biológica.

### Resultado de conjuntos gênicos

Houve conjuntos FDR-significativos, incluindo termos do GWAS Catalog para insônia, sintomas de insônia e curta duração do sono. Esses resultados fornecem coerência fenotípica para a lista, mas não são validação totalmente independente, pois catálogos de GWAS podem incluir estudos relacionados.

Os termos moleculares de ligação ao RNA 7SK e inibição de quinases dependentes de ciclina foram conduzidos por HEXIM1 e HEXIM2. A interpretação deve mencionar explicitamente que esses termos dependem de dois genes. O conjunto posicional 17q21 também foi significativo e deve ser relatado como concentração em locus, não como via biológica independente.

## 9. Arquivos principais gerados

- `Resultados/analise_padronizada/11_provisional_consensus_protein_coding.tsv`: lista dos 16 genes e suas evidências.
- `Resultados/analise_padronizada/20_core16_x_opengwas_ukb_magma.tsv`: resultados gene a gene nos dois UKB de sensibilidade.
- `Resultados/analise_padronizada/21_core16_opengwas_ukb_sensitivity_summary.tsv`: resumo da sensibilidade UKB.
- `Resultados/analise_padronizada/22_core16_final_evidence_table.tsv`: tabela final integrada para as próximas etapas.
- `Resultados/analise_padronizada/23_gene2func_fdr_significant_gene_sets.tsv`: conjuntos significativos no GENE2FUNC.
- `Resultados_GENE2FUNC/FUMA_gene2func761816_16_genes/`: resultados completos do GENE2FUNC.

## 10. Limitações a declarar

- Fenótipos de insônia diferem entre estudos; Lane representa fenótipo extremo, enquanto outros estudos representam definições mais amplas.
- A indisponibilidade dos dados 23andMe reduz o tamanho efetivo das análises públicas de Jansen e Watanabe.
- Jansen, Watanabe e os dois conjuntos UKB OpenGWAS não devem ser tratados como quatro replicações independentes.
- Mapeamento SNP2GENE prioriza genes, mas não identifica causalidade molecular.
- Resultados de enrichment de uma lista pequena podem ser influenciados por poucos genes, por blocos de LD e por proximidade genômica.

## 11. Próxima etapa prevista

Usar a tabela integrada de 16 genes para análises complementares em g:Profiler, ShinyGO, KEGG/Reactome e STRING. Essas ferramentas devem ser apresentadas como interpretação funcional posterior da mesma lista priorizada, sem redefinir o núcleo de genes com base em termos de vias obtidos depois.
