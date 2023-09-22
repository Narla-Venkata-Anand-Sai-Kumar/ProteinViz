from flask import Flask, render_template, request, redirect, url_for
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

app = Flask(__name__)

# Define the route for the homepage
@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        smiles = request.form['smiles']
        mol = Chem.MolFromSmiles(smiles)

        if mol is not None:
            AllChem.EmbedMolecule(mol, randomSeed=42)

            viewer = py3Dmol.view(width=500, height=500)
            viewer.addModel(Chem.MolToMolBlock(mol), "mol")
            viewer.setStyle({"stick": {}})
            viewer.zoomTo()

            mol3d_html = viewer.render()
            html, u = viewer._make_html()
            print(u)
            html = f'''
            <!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>3D Protein Visualization</title>
    <!-- Link to Bootstrap CSS -->
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.5.2/css/bootstrap.min.css">
</head>
<body style="margin: 0; padding: 0; font-family: Arial, sans-serif; background: linear-gradient(to bottom, #9900cc, #6600cc, #0033cc); color: #fff; line-height: 1.6; text-align: center; background-attachment: fixed;">

    <!-- Navigation bar using Bootstrap -->
    <nav class="navbar navbar-expand-lg navbar-light bg-light" style="background-color: transparent !important;">
        <a class="navbar-brand" href="#" style="font-size: 24px; color: #fff !important;">ProteinViz</a>
        <button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarNav" aria-controls="navbarNav" aria-expanded="false" aria-label="Toggle navigation">
            <span class="navbar-toggler-icon" style="background-color: #fff !important;"></span>
        </button>
        <div class="collapse navbar-collapse" id="navbarNav">
            <ul class="navbar-nav ml-auto">
                <li class="nav-item"><a class="nav-link" href="#" style="color: #fff !important;">Home</a></li>
                <li class="nav-item"><a class="nav-link" href="#" style="color: #fff !important;">About</a></li>
                <li class="nav-item"><a class="nav-link" href="#" style="color: #fff !important;">Services</a></li>
                <li class="nav-item"><a class="nav-link" href="#" style="color: #fff !important;">Contact</a></li>
            </ul>
        </div>
    </nav>

    <!-- Hero section using Bootstrap -->
    <section class="hero text-center" style="background-image: url('background-image.jpg'); background-size: cover; background-position: center; padding: 100px 0; text-shadow: 2px 2px 4px rgba(0, 0, 0, 0.5);">
        <div class="container">
            <h1 class="display-4" style="font-size: 36px;">Welcome to Our 3D Protein Visualization Page</h1>
            <p class="lead" style="font-size: 24px;">Uncover the beauty of molecular structures in 3D</p>
        </div>
        <h1 style="font-size: 36px;">3D visualization of {smiles}</h1>
    </section>
    <section id="3dview" style="display: flex; justify-content: center; align-items: center; height: 60vh;">
        {html}
    </section>

    <footer class="bg-light py-4" style="background-color: transparent !important; color: #fff; padding: 10px 0;">
        <div class="container text-center">
            <div class="footer-logo" style="font-size: 18px;">
                <p>Developed By <b>N V ANAND SAI KUMAR</b></p>
            </div>
            <p>&copy; 2023 3D Protein Visualization</p>
            <p><a href="#">Privacy Policy</a> | <a href="#">Terms of Service</a></p>
        </div>
    </footer>

    <!-- Link to Bootstrap JavaScript Scripts -->
    <script src="https://code.jquery.com/jquery-3.5.1.slim.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/@popperjs/core@2.5.3/dist/umd/popper.min.js"></script>
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.5.2/js/bootstrap.min.js"></script>
</body>
</html>
'''
            f = open('templates/main.html',"wt")
            f.write(html)
            # Redirect to the output page with the visualization
            return redirect(url_for('output', smiles=smiles, mol3d_html=mol3d_html))
        else:
            error_message = "Invalid SMILES input. Please try again."
            return render_template('index.html', error_message=error_message)

    return render_template('index.html')

# Define the route for the output page
@app.route('/output')
def output():
    smiles = request.args.get('smiles')
    mol3d_html = request.args.get('mol3d_html')
    return render_template('main.html')

if __name__ == '__main__':
    app.run(debug=True)
