# Rede STRING: lista principal e parâmetros

## Arquivo de entrada

Use `core16_string_symbols.txt`. Ele contém somente os 16 genes núcleo protein-coding, um símbolo HGNC por linha.

## Submissão

1. Acesse `https://string-db.org/`.
2. Selecione **Multiple proteins**.
3. Cole o conteúdo do arquivo de entrada.
4. Em organismo, escolha *Homo sapiens*.
5. Confirme o mapeamento de todos os 16 genes antes de continuar.

## Configuração principal

- máximo de interatores adicionais: 0;
- network type: full STRING network;
- network flavor: confidence;
- confidence score: 0,70;
- mantenha todos os canais de evidência ativos;
- não remova nós desconectados; eles são resultado biologicamente informativo.

Para a síntese de termos funcionais no relatório, foram preservados apenas termos com FDR <= 0,01. Esse filtro de reporte não altera a rede nem o critério de confiança das arestas.

## O que salvar

- imagem em SVG e PNG;
- tabela de interações, com combined score e canais de evidência;
- resumo de network statistics, especialmente número de nós, arestas, grau médio, coeficiente de clustering e PPI enrichment p-value;
- tabela de functional enrichment, mesmo que não haja termos significativos;
- captura das configurações e da versão do STRING.

## Interpretação

O STRING integra evidências de associação funcional/proteica. Uma aresta não prova interação física direta nem causalidade para insônia. A rede não deve receber proteínas adicionais na figura principal, pois esses nós não foram selecionados pelo pipeline GWAS.
