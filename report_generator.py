from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from datetime import datetime

def generate_pdf_report(result, filename="patent_report.pdf"):
    """
    Generates a detailed PDF report based on the selected patent result.
    """

    doc = SimpleDocTemplate(filename, pagesize=letter)
    styles = getSampleStyleSheet()
    story = []

    title = f"<b>Patent Report</b>"
    story.append(Paragraph(title, styles['Title']))
    story.append(Spacer(1, 20))

    fields = {
        "SMILES": result.get("smiles", ""),
        "Patent ID": result.get("patent_id", ""),
        "Title": result.get("title", ""),
        "Abstract": result.get("abstract", ""),
        "URL": result.get("url", ""),
        "Publication Year": result.get("year", "")
    }

    for key, val in fields.items():
        text = f"<b>{key}:</b> {val}"
        story.append(Paragraph(text, styles['BodyText']))
        story.append(Spacer(1, 12))

    timestamp = f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
    story.append(Paragraph(timestamp, styles['Italic']))

    doc.build(story)

    return filename
