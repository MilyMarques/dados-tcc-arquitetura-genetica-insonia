# Metodologia pós-GENE2FUNC e plano da rede STRING

## Escopo

> **Nota de atualização.** Este documento começou como plano de análise. Os parâmetros e resultados finais do STRING são os registrados em `DOCUMENTACAO_STRING_E_INTERPRETACAO.md`: escore combinado >=0,700, sem interatores adicionais e termos reportados com FDR <=0,01. Onde houver divergência, prevalece a documentação final.

Este documento registra as análises funcionais posteriores à definição do núcleo de 16 genes. Ele não altera a regra de seleção do núcleo, que já foi documentada em `DOCUMENTACAO_METODOLOGIA_ATE_GENE2FUNC.md`.

## 1. Princípio analítico

Foram aplicadas duas famílias de análise funcional, que não devem ser confundidas:

1. **Sobre-representação de uma lista de genes:** pergunta se os genes selecionados ocorrem em um termo/via mais frequentemente que o esperado no universo de background.
2. **MAGMA gene-set analysis competitivo:** pergunta se genes pertencentes a um conjunto pré-definido apresentam, como grupo, associações gene-based mais fortes com insônia do que genes fora daquele conjunto.

A primeira família foi aplicada ao núcleo de 16 genes e, como sensibilidade histórica, à união de 28 genes MAGMA de Jansen/Watanabe. A segunda usa o sinal de todo o genoma de cada GWAS e é independente da lista core.

## 2. GENE2FUNC

### Entrada e teste

- Lista principal: 16 genes núcleo, identificadores Ensembl.
- Background: 18.717 genes MAGMA-testáveis nos quatro estudos principais.
- Teste: hipergeométrico de sobre-representação.
- Correção: Benjamini-Hochberg, FDR < 0,05, realizada pelo FUMA.
- Bases de expressão: GTEx v8 (53 tecidos e 30 categorias gerais) e BrainSpan.
- MHC: não reaplicado no GENE2FUNC, pois o teste recebe genes e background, não loci.

### Interpretação

Não houve enriquecimento de expressão diferencial por tecido após FDR. Houve conjuntos gênicos FDR-significativos no FUMA, incluindo termos do GWAS Catalog relacionados a insônia; estes foram tratados como coerência fenotípica e não como validação independente.

## 3. g:Profiler g:GOSt

### Entrada e universo

- Organismo: *Homo sapiens*.
- Lista principal: 16 genes Ensembl.
- Lista de sensibilidade: 28 genes MAGMA históricos de Jansen/Watanabe.
- Background principal: 18.717 genes MAGMA-testáveis.
- Background de sensibilidade: 19.141 genes com p MAGMA finito nas duas rodagens preliminares.
- Tipo de lista: não ordenada.
- Fontes primárias: GO:BP, GO:MF, GO:CC, KEGG, Reactome e WikiPathways.
- Domínio estatístico: `Custom over annotated genes`.

### Estatística

Para cada termo, o g:GOSt calcula um teste hipergeométrico unilateral/Fisher para sobre-representação. O p-valor bruto é a probabilidade de observar a sobreposição encontrada, ou maior, se genes fossem sorteados ao acaso do background.

Os p-valores são corrigidos dentro da análise g:Profiler pelo **g:SCS** a 0,05. g:SCS é o procedimento padrão do g:Profiler e foi criado para acomodar a estrutura dependente e hierárquica de termos funcionais, especialmente do Gene Ontology. Essa correção é independente das correções realizadas anteriormente no MAGMA, porque corrige outra família de testes.

`Custom over annotated genes` significa que, para cada fonte, o denominador contém somente genes do background customizado que possuem pelo menos uma relação gene-termo naquela fonte. Por exemplo, no GO:MF, genes sem anotação GO:MF não entram no domínio daquele teste. As anotações são relações curadas/integradas entre genes e termos em GO, KEGG, Reactome e WikiPathways.

### Resultado principal e sensibilidade

Não houve termo significativo após g:SCS para a lista núcleo de 16 genes ou para a lista histórica de 28 genes. Resultados mostrados por “All results” são nominais e não são interpretados como achados confirmatórios.

## 4. ShinyGO

### Entrada e estatística

- Organismo: humano.
- Mesmas listas e respectivos backgrounds do g:Profiler.
- Limiar: FDR Benjamini-Hochberg < 0,05.
- Não foi selecionada a opção “Use pathway DB for gene counts”, preservando o background MAGMA customizado.
- Tamanho mínimo de via: 2 genes; máximo: 5.000 genes.
- Remoção de redundância ativada para apresentação de termos GO.

O ShinyGO aplica teste hipergeométrico e corrige os p-valores por Benjamini-Hochberg. O FDR do ShinyGO não é diretamente equivalente ao g:SCS do g:Profiler; as plataformas foram usadas para verificar coerência temática, não para combinar p-valores.

### Resultado

Não houve enrichment significativo após FDR para a lista núcleo de 16 genes nem para a sensibilidade histórica de 28 genes. A identificação dos genes de entrada foi conferida no arquivo `geneInfo` exportado em cada rodada.

## 5. MAGMA competitive gene-set analysis

### Entrada e estatística

Foram extraídos os arquivos `magma.gsa.out` dos quatro jobs FUMA principais: Hammerschlag 2017, Lane 2019, Jansen 2019 público e Watanabe 2022 público. Cada arquivo continha 17.010 conjuntos gênicos.

O MAGMA aplica uma regressão competitiva gene-set, testando se a associação gene-based média de genes do conjunto é maior que a de genes fora do conjunto, enquanto controla covariáveis do MAGMA, incluindo tamanho e densidade gênica e medidas relacionadas ao número de SNPs. O teste é unilateral, na direção positiva do efeito do conjunto.

Para cada estudo foram recalculados:

- FDR Benjamini-Hochberg sobre todos os 17.010 conjuntos;
- Bonferroni: 0,05 / 17.010 = 2,94 x 10^-6.

Jansen e Watanabe foram considerados uma única unidade `ukb_public` na avaliação de recorrência entre estudos, devido à sobreposição de participantes UK Biobank.

### Regra de interpretação

- FDR < 0,05 em ao menos uma unidade: resultado reportável, porém exploratório se presente em apenas uma unidade.
- FDR < 0,05 em pelo menos duas unidades analíticas: resultado consistente entre análises; a sobreposição parcial de participantes deve ser considerada na interpretação.

### Resultado

Hammerschlag, Jansen e Watanabe não tiveram gene sets FDR-significativos. Lane teve três conjuntos FDR-significativos, sem significância Bonferroni e sem recorrência em outra unidade. Portanto, não foram identificadas vias MAGMA consistentemente associadas à insônia em duas ou mais unidades analíticas.

## 6. Comparação externa com Watanabe (2022)

As Tabelas Suplementares 14 e 37 são usadas como benchmark de literatura:

- Tabela 14: genes mapeados por FUMA no estudo publicado;
- Tabela 37: genes classificados pelos autores como alta confiança, incluindo HCP.

Essa comparação não constitui replicação independente, pois Watanabe usa componente UK Biobank sobreposto, mas documenta concordância metodológica e bibliográfica.

## 7. STRING: configuração final

A rede principal conterá somente os 16 genes núcleo protein-coding. A lista será enviada como símbolos HGNC para *Homo sapiens*.

Parâmetros a registrar no momento da execução:

- versão do STRING exibida na interface;
- espécie: *Homo sapiens*;
- máximo de interatores adicionais: 0;
- rede: full STRING network, com todos os canais de evidência ativos;
- visualização: confidence;
- limiar de confiança: 0,70 (alto);
- análise de enriquecimento PPI fornecida pelo STRING;
- exportações: SVG/PNG, tabela de interações e tabela de enrichment funcional.

O score combinado do STRING agrega evidências de experimentos, bases curadas, coexpressão, vizinhança gênica, fusão, coocorrência e text-mining. Ele indica confiança de associação funcional/proteica, não prova interação física direta nem causalidade em insônia.

Os termos funcionais resumidos na tabela integrada foram restringidos a FDR <=0,01. Não foram adicionados nós externos na figura principal.

## 8. Conclusão metodológica provisória

O núcleo de 16 genes é sustentado por associação gene-based MAGMA e mapeamento SNP2GENE, com consistência adicional nos dois estudos UKB de sensibilidade. As análises posteriores de sobre-representação e gene sets foram interpretativas e não modificaram o núcleo. A ausência de vias replicadas após correção múltipla é um resultado válido, compatível com arquitetura poligênica e heterogeneidade fenotípica, sem demonstrar ausência absoluta de mecanismos biológicos compartilhados.
