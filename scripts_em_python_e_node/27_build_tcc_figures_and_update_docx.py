"""Create the final Core16/STRING figures and insert them into a copy of the TCC DOCX.

The R counterpart for the Core16 heatmap is R/21_plot_core16_evidence_heatmap.R.
This script creates the same evidence matrix with Pillow so the document can be rendered
in this workspace, where an R runtime is unavailable.
"""
from pathlib import Path
import csv
import shutil
from PIL import Image, ImageChops, ImageDraw, ImageFont
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt
from docx.text.paragraph import Paragraph
from docx.oxml import OxmlElement

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT / "projeto_insonia_tcc"
FIGDIR = PROJECT / "Figuras_TCC"
ANALYTIC = PROJECT / "Resultados" / "analise_padronizada"
SOURCE = ROOT / "TCC_Arquitetura_Genetica_Insônia_Discussao.docx"
OUTPUT = ROOT / "TCC_Arquitetura_Genetica_Insônia_versao_figuras_revisada.docx"

def font(size, bold=False):
    candidates = [
        Path("C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf"),
        Path("C:/Windows/Fonts/calibrib.ttf" if bold else "C:/Windows/Fonts/calibri.ttf"),
    ]
    for p in candidates:
        if p.exists():
            return ImageFont.truetype(str(p), size)
    return ImageFont.load_default()

def crop_white(image, padding=25):
    rgb = image.convert("RGB")
    bg = Image.new("RGB", rgb.size, "white")
    diff = ImageChops.difference(rgb, bg)
    bbox = diff.getbbox()
    if not bbox:
        return rgb
    left, top, right, bottom = bbox
    return rgb.crop((max(0, left-padding), max(0, top-padding), min(rgb.width, right+padding), min(rgb.height, bottom+padding)))

def contain(im, size):
    im = im.copy()
    im.thumbnail(size, Image.Resampling.LANCZOS)
    return im

def create_panel():
    core = crop_white(Image.open(PROJECT / "Resultados_STRING" / "STRING_core16" / "STRING_network_core16.png"))
    legacy = crop_white(Image.open(PROJECT / "Resultados_STRING" / "STRING_legacy28" / "STRING_network_legacy28.png"))
    global_net = crop_white(Image.open(PROJECT / "Resultados_STRING" / "STRING_global_gwas_catalog" / "STRING_network_global_gwas_catalog.png"))
    canvas = Image.new("RGB", (3600, 4050), "white")
    draw = ImageDraw.Draw(canvas)
    title_font, header_font, note_font = font(72, True), font(42, True), font(28, False)
    draw.text((1800, 80), "Redes STRING em três níveis de evidência", fill="#000000", font=title_font, anchor="ma")
    draw.text((1800, 180), "Homo sapiens | full STRING network | score ≥ 0,700 | sem nós adicionais", fill="#404040", font=note_font, anchor="ma")
    # Panel A: deliberately larger because it represents the broad descriptive layer.
    draw.rounded_rectangle((90, 260, 3510, 2180), radius=20, outline="#B7C9D6", width=4, fill="#F9FBFC")
    draw.text((150, 305), "A", fill="#1F4E79", font=header_font)
    draw.text((230, 315), "GWAS Catalog global — camada descritiva", fill="#1F1F1F", font=header_font)
    g = contain(global_net, (3250, 1720))
    canvas.paste(g, ((3600-g.width)//2, 430))
    # Panels B and C: the sensitivity and the stringent primary core.
    for x, letter, label, im in [(90, "B", "Legacy28 — sensibilidade histórica", legacy), (1830, "C", "Core16 — conjunto principal", core)]:
        draw.rounded_rectangle((x, 2260, x+1680, 3800), radius=20, outline="#B7C9D6", width=4, fill="#F9FBFC")
        draw.text((x+60, 2310), letter, fill="#1F4E79", font=header_font)
        draw.text((x+140, 2320), label, fill="#1F1F1F", font=header_font)
        fitted = contain(im, (1540, 1280))
        canvas.paste(fitted, (x+(1680-fitted.width)//2, 2460))
    draw.text((1800, 3915), "A conectividade cresce com a amplitude da lista; a rede global não redefine o Core16.", fill="#404040", font=note_font, anchor="ma")
    output = FIGDIR / "painel_redes_STRING.png"
    canvas.save(output, dpi=(300, 300), optimize=True)
    return output

def create_heatmap():
    with open(ANALYTIC / "22_core16_final_evidence_table.tsv", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    cols = ["MAGMA\nunidades", "Posicional", "eQTL\ncerebral", "Cromatina\ncortical"]
    width, height, left, top, cell_w, cell_h = 2200, 2000, 520, 350, 360, 80
    im = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(im)
    title_font, label_font, cell_font, note_font = font(52, True), font(32, True), font(28, False), font(24, False)
    draw.text((width//2, 70), "Evidências integradas que definem o Core16", font=title_font, fill="#000000", anchor="ma")
    draw.text((width//2, 145), "MAGMA: número de unidades independentes com significância por Bonferroni", font=note_font, fill="#404040", anchor="ma")
    for j, col in enumerate(cols):
        cx = left + j*cell_w + cell_w//2
        draw.multiline_text((cx, top-110), col, font=label_font, fill="#1F1F1F", anchor="ma", align="center", spacing=3)
    for i, r in enumerate(rows):
        y = top + i*cell_h
        draw.text((left-25, y+cell_h//2), r["SYMBOL"], font=label_font, fill="#1F1F1F", anchor="rm")
        values = [int(r["n_independent_units_magma_bonf"]), r["positional"] == "TRUE", r["eqtl_brain"] == "TRUE", r["chromatin_cortex"] == "TRUE"]
        for j, value in enumerate(values):
            if j == 0:
                color = "#1F4E79" if value == 3 else "#5B9BD5"
                text = str(value)
            else:
                color = "#5B9BD5" if value else "#F1F1F1"
                text = "Sim" if value else "Não"
            x = left + j*cell_w
            draw.rectangle((x, y, x+cell_w-3, y+cell_h-3), fill=color, outline="white", width=2)
            draw.text((x+cell_w//2, y+cell_h//2), text, font=cell_font, fill="white" if color == "#1F4E79" else "black", anchor="mm")
    draw.text((width//2, top + len(rows)*cell_h + 70), "Colunas de mapeamento representam modalidades SNP2GENE presentes em pelo menos uma unidade analítica da síntese.", font=note_font, fill="#404040", anchor="ma")
    output = FIGDIR / "core16_evidence_heatmap.png"
    im.save(output, dpi=(300,300), optimize=True)
    return output

def insert_after(paragraph, style=None):
    new_p = OxmlElement("w:p")
    paragraph._p.addnext(new_p)
    return Paragraph(new_p, paragraph._parent)

def delete_paragraph(paragraph):
    parent = paragraph._p.getparent()
    parent.remove(paragraph._p)

def set_caption(paragraph, text):
    paragraph.style = "Normal"
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    paragraph.paragraph_format.space_before = Pt(3)
    paragraph.paragraph_format.space_after = Pt(9)
    run = paragraph.add_run(text)
    run.font.size = Pt(9)

def add_figure_after(paragraph, image_path, width_inches, caption):
    fig_p = insert_after(paragraph)
    fig_p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    fig_p.paragraph_format.space_before = Pt(9)
    fig_p.paragraph_format.space_after = Pt(3)
    fig_p.add_run().add_picture(str(image_path), width=Inches(width_inches))
    cap = insert_after(fig_p)
    set_caption(cap, caption)
    return cap

def table_text(table):
    return " ".join(cell.text for row in table.rows for cell in row.cells)

def revise_table4(doc):
    table = next(t for t in doc.tables if "Core16 com FDR" in table_text(t))
    header = ["Análise", "Genes MAGMA\navaliados", "Core16 com FDR", "Core16 com Bonferroni", "Interpretação"]
    for cell, value in zip(table.rows[0].cells, header): cell.text = value
    values = [
        ["OpenGWAS UKB - execução 1", "18.774", "16", "13", "Sensibilidade interna"],
        ["OpenGWAS UKB - execução 2", "18.838", "16", "16", "Sensibilidade interna"],
        ["FinnGen R13 corrigido", "18.802", "Não avaliado", "0", "Validação externa conservadora"],
    ]
    for row, vals in zip(table.rows[1:], values):
        for cell, value in zip(row.cells, vals): cell.text = value

def remove_separate_string_figures(doc):
    targets = ["Rede de associação funcional do Core16", "Rede de associação funcional da lista Legacy28", "Rede de associação funcional do catálogo global"]
    body = doc._body._element
    for paragraph in list(doc.paragraphs):
        if any(paragraph.text.startswith(t) for t in targets):
            prev = paragraph._p.getprevious()
            if prev is not None and (prev.xpath('.//w:drawing') or prev.xpath('.//w:pict')):
                body.remove(prev)
            delete_paragraph(paragraph)

def find_paragraph(doc, starts):
    for p in doc.paragraphs:
        if p.text.startswith(starts): return p
    raise ValueError(f"Paragraph not found: {starts}")

def update_document(panel, heatmap, graph):
    shutil.copy2(SOURCE, OUTPUT)
    doc = Document(OUTPUT)
    revise_table4(doc)
    remove_separate_string_figures(doc)
    core_source = find_paragraph(doc, 'Fonte: elaboração própria. “Sim” indica')
    core_caption = "Mapa de calor das evidências integradas que definem o Core16. A coluna MAGMA informa o número de unidades independentes com significância por Bonferroni; as demais colunas mostram modalidades de mapeamento SNP2GENE. Fonte: elaboração própria."
    add_figure_after(core_source, heatmap, 6.2, core_caption)
    string_source = find_paragraph(doc, 'Fonte: elaboração própria a partir do STRING. As arestas')
    panel_caption = "Painel comparativo das redes de associação funcional no STRING: (A) catálogo global do GWAS Catalog, camada descritiva; (B) Legacy28, análise de sensibilidade histórica; e (C) Core16, conjunto principal priorizado. Fonte: elaboração própria a partir do STRING."
    panel_cap = add_figure_after(string_source, panel, 6.35, panel_caption)
    graph_caption = "Grafo de processos funcionais enriquecidos na rede global de genes associados à insônia no STRING. Os agrupamentos representam termos de processos anotados; a análise é descritiva e não demonstra causalidade. Fonte: elaboração própria a partir do STRING."
    add_figure_after(panel_cap, graph, 6.1, graph_caption)
    doc.save(OUTPUT)
    return OUTPUT

def main():
    FIGDIR.mkdir(exist_ok=True)
    heatmap = create_heatmap()
    panel = create_panel()
    graph_source = PROJECT / "Resultados_STRING" / "STRING_global_gwas_catalog" / "enrichment_Process_sim0.8_graph_global_gwas_catalog.png"
    graph = FIGDIR / "grafo_processos_STRING_global.png"
    shutil.copy2(graph_source, graph)
    output = update_document(panel, heatmap, graph)
    print(f"Created: {heatmap}")
    print(f"Created: {panel}")
    print(f"Created: {graph}")
    print(f"Created: {output}")

if __name__ == "__main__":
    main()
