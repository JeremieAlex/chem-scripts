import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageDraw
from fpdf import FPDF

import os
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageDraw
from fpdf import FPDF

def save_smiles_to_pdf(df, smiles_column, name_column, output_pdf):
    # Sort the DataFrame based on 'norm_value'

    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=15)

    structures_per_page = 1
    page_width = 600  # Width of A3 page in mm (297 mm x 420 mm)
    page_height = 600  # Height of A3 page in mm (420 mm x 594 mm)
    structure_width = (page_width)/3.5  # Divide page width into 2 for 2 structures per row
    structure_height = (page_height)/3.5  # Divide page height into 2 for 2 rows

    for i, (smiles, name) in enumerate(zip(df[smiles_column], df[name_column])):
        # Generate RDKit molecule object from SMILES
        mol = Chem.MolFromSmiles(smiles)

        if mol:
            # Generate a 2D image of the molecule
            img = Draw.MolToImage(mol, size=(400, 400))

            # Convert PIL image to grayscale
            img_gray = img.convert('RGB')

            # Create new image with grayscale data
            img_new = Image.new('RGB', img_gray.size)
            img_new.paste(img_gray)

            # Create a drawing context
            draw = ImageDraw.Draw(img_new)

            # Add compound name to the image
            draw.text((10, 10), name, fill='black')  # Adjust position as needed

            # Save the image as a PNG file
            img_path = f"{i}.png"
            img_new.save(img_path)

            # Add the PNG image to the PDF
            if i % structures_per_page == 0:
                pdf.add_page()
                x = 10
                y = 10

            pdf.image(img_path, x=x, y=y, w=structure_width, h=structure_height)

            # Update x and y coordinates for next structure
            x += structure_width + 10
            if (i + 1) % 2 == 0:
                x = 10
                y += structure_height + 10

            # Cleanup: Remove the temporary image file
            os.remove(img_path)

    # Save the PDF
    pdf.output(output_pdf)

# Example usage:
df = pd.read_csv('12hr_screens_Feb2024/2sd_hit_compounds.csv')
#df = df[['Compound Name', 'SMILES']]
# Convert the SMILES column to strings
df['SMILES'] = df['SMILES'].astype(str)
df.sort_values(by='mean_norm', inplace=True, ascending=True)


# Example usage:
output_pdf = "2SD_hits_molecules_ordered.pdf"
save_smiles_to_pdf(df, 'SMILES', 'Compound Name', output_pdf)
