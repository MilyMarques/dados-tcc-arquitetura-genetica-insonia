# Guia prático: g:Profiler e ShinyGO

## Objetivo e regra de consistência

Estas duas ferramentas serão usadas para interpretar funcionalmente o mesmo conjunto principal de 16 genes. Elas **não** devem ser usadas para adicionar/remover genes da lista núcleo. Os resultados são complementares ao GENE2FUNC porque usam bases, versões e métodos próprios.

Use sempre:

- lista de interesse: `enrichment_inputs/core16_ensembl_ids.txt`;
- background: `enrichment_inputs/magma_background_18717_ensembl_ids.txt`;
- espécie: *Homo sapiens*;
- análise não ordenada (lista simples, sem ranking);
- foco interpretativo: GO, Reactome, KEGG e WikiPathways.

Use os identificadores Ensembl como entrada primária: assim, o background mantém exatamente os 18.717 genes testáveis pelo MAGMA. Os arquivos com símbolos HGNC são auxiliares para leitura humana e para a etapa posterior no STRING; alguns genes testáveis não possuem símbolo disponível, por isso essa versão tem menos linhas.

## 1. g:Profiler (g:GOSt)

Acesse `https://biit.cs.ut.ee/gprofiler/gost`.

1. Em **Organism**, selecione *Homo sapiens*.
2. Abra `core16_ensembl_ids.txt`, copie as 16 linhas e cole em **Query**.
3. Em **Statistical domain scope**, selecione **Custom over annotated genes**.
4. Cole ou carregue as linhas de `magma_background_18717_ensembl_ids.txt` no campo de **custom background**.
5. Mantenha a query como não ordenada. Não use “ordered query”.
6. Em fontes de dados, marque somente:
   - GO: Biological Process (GO:BP);
   - GO: Molecular Function (GO:MF);
   - GO: Cellular Component (GO:CC);
   - Reactome (REAC);
   - KEGG;
   - WikiPathways (WP).
7. Não marque GWAS Catalog, HPO, motivos de fator de transcrição, miRNA, CORUM ou Human Protein Atlas nesta rodada primária. Isso reduz testes exploratórios e evita repetir diretamente os achados de GWAS de insônia já vistos no FUMA.
8. Mantenha **g:SCS** como correção múltipla e limiar 0,05. g:SCS é a correção padrão da ferramenta e lida com a dependência entre termos GO. Não misture esse resultado com FDR como se fossem o mesmo teste.
9. Mantenha IEA incluído nesta primeira rodada; registrar essa decisão é suficiente.
10. Execute a análise. Baixe a tabela completa e a figura Manhattan plot/GO context, se houver resultados significativos.

### Como registrar os resultados

Salve a exportação na pasta `Resultados_gProfiler` com um nome que inclua a data. Para cada termo a relatar, registre fonte, ID, nome do termo, p ajustado, tamanho do termo, número de genes sobrepostos e quais genes da lista explicam o sinal.

Se nenhum resultado permanecer significativo, este é um resultado válido: a lista pequena não apresenta enriquecimento detectável nessas bases e sob esse background. Não altere o limiar de significância para produzir resultados.

## 2. ShinyGO

Acesse `https://bioinformatics.sdstate.edu/go/`.

1. Selecione *Homo sapiens*.
2. Cole as 16 linhas de `core16_ensembl_ids.txt` no campo da lista de genes.
3. Abra a opção **Customized background genes** e carregue/cole `magma_background_18717_ensembl_ids.txt`.
4. Não selecione **Use pathway database for gene counts**. O background personalizado MAGMA deve ser mantido como universo estatístico.
5. Estabeleça FDR < 0,05 e mantenha a correção Benjamini-Hochberg da ferramenta.
6. Selecione as bases GO:BP, GO:MF, GO:CC, Reactome, KEGG e WikiPathways, se apresentadas na interface atual.
7. Para a visualização, habilite **Remove redundancy** em resultados GO. Isso reduz termos quase idênticos; guarde também a tabela completa antes de resumir.
8. Em “sort pathway by”, comece com **FDR**. Depois examine Fold Enrichment como medida de tamanho de efeito; não escolha termos somente pelo maior fold enrichment se forem sustentados por apenas dois genes.
9. Baixe: tabela de enrichment, gráfico, árvore de termos e rede de termos, se existirem resultados significativos.

Não use ainda a aba STRING do ShinyGO. A rede final deve ser feita diretamente no STRING, em uma etapa posterior, para preservar uma configuração documentável e imagens de melhor qualidade.

## 3. Como interpretar e comparar

- Resultado significativo: p ajustado < 0,05 segundo o método daquela plataforma.
- Priorize termos de vias/GO sustentados por pelo menos três genes da lista, sempre que possível.
- Verifique os genes que sustentam cada termo. `HEXIM1` e `HEXIM2`, por exemplo, podem gerar termos em conjunto; isso deve ser explicitado, não generalizado como uma via ampla.
- Termos muito parecidos devem ser agrupados em um tema biológico, em vez de relatados como descobertas independentes.
- Concordância entre FUMA GENE2FUNC, g:Profiler e ShinyGO é coerência funcional; não é validação independente, pois as ferramentas usam bancos de dados parcialmente sobrepostos.
- Diferenças entre plataformas são esperadas, por causa de bases de dados, versões, conversão de IDs e correção múltipla diferentes.

## 4. O que trazer para a próxima interpretação

Depois de cada plataforma, coloque a pasta/arquivo exportado dentro do projeto. A integração posterior deve reter somente termos que atendam ao limiar ajustado e registrar, para cada um, os genes sobrepostos, banco de origem e método de correção.
