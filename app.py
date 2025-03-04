from flask import Flask, render_template, request, redirect, url_for
from rdkit import Chem
from rdkit.Chem import AllChem
import require as py3Dmol

app = Flask(__name__)
def codon_category(codon):
     codon_dict={"UUU":"Phe","UUC":"Phe","UUA":"Leu","UUA":"Leu",
                 "UCU":"Ser","UCC":"Ser","UCA":"Ser","UCG":"Ser",
                 "UAU":"Tyr","UAC":"Tyr","UAA":"Stop","UAG":"Stop",
                 "UGU":"Cys","UGC":"Cys","UGA":"Stop","UGG":"Trp",
                 "CUU":"Leu","CUC":"Leu","CUA":"Leu","CUG":"Leu",
                 "CCU":"Pro","CAC":"Pro","CCA":"Pro","CCG":"Pro",
                 "CAU":"His","CAC":"His","CAA":"Gln","CAG":"Gln",
                 "CGU":"Arg","CGC":"Arg","CGA":"Arg","CGG":"Arg",
                 "AUU":"Ile","AUC":"Ile","AUA":"Ile","AUG":"Met",
                 "ACU":"Thr","ACC":"Thr","ACA":"Thr","ACG":"Thr",
                 "AAU":"Asn","AAC":"Asn","AAA":"Lys","AAG":"Lys",
                 "AGU":"Ser","AGC":"Ser","AGA":"Arg","AGG":"Arg",
                 "GUU":"Val","GUC":"Val","GUA":"Val","GUG":"Val",
                 "GCU":"Ala","GCC":"Ala","GCA":"Ala","GCA":"Ala",
                 "GAU":"Asp","GAC":"Asp","GAA":"Glu","GAG":"Glu",
                 "GGU":"Gly","GGC":"Gly","GGA":"Gly","GGG":"Gly"}
     return codon_dict[codon]
def amino_acid(sequence):
     sequence_length=len(sequence)-3
     amino=""
     for i in range(0,sequence_length,3):
         code=str(sequence[i]+sequence[i+1]+sequence[i+2])
         amino+=codon_category(code)
         amino+='-'
     return amino[:-1]
def Smile_String(smile):
     Smile_dict={"Ala":"CC(C(=O)O)N","Cys":" CSCC(C(=O)O)N","Asp":"CC(C(=O)O)C(=O)O",
                 "Glu":"CC(C(=O)O)C(=O)OC(C)C","Phe":"Nc1ccccc1C(=O)O",
                 "Gly":"NCC(=O)O","His":"Nc1c[nH]cn1C(=O)O","Ile":"CCC(C)C(C(=O)O)N",
                 "Lys":"NCCCCC@HN","Leu":"CCC(C)CC(C(=O)O)N","Met":"CSCCC(C(=O)O)N",
                 "Asn":"CC(C(=O)O)C(=O)N","Pro":"CC1CC(C(=O)N)C1","Gln":"CC(C(=O)O)C(=O)NC",
                 "Arg":"N=C(N)NCCCCC@HN","Ser":"OCC(C(=O)O)N","Thr":"CC(C(=O)O)C(O)N",
                 "Val":"CCC(C)C(C(=O)O)N","Trp":"Nc1ccc2c(c1)CC@@HN2","Tyr":"Nc1ccc(O)cc1C(=O)O"}
     return Smile_dict[smile]
 
# Define the route for the homepage
@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        smiles = request.form['smiles']
        dump_smile = smiles
        k=amino_acid(smiles)
        S=k.split('-')
        print(S)
        smiles = ""
        for i in S:
            smiles+=Smile_String(i)
        print(smiles)
        mol = Chem.MolFromSmiles(smiles)

        if mol is not None:
            
            AllChem.EmbedMolecule(mol, randomSeed=42)

            viewer = py3Dmol.view(width=500, height=500)
            viewer.addModel(Chem.MolToMolBlock(mol), "mol")
            viewer.setStyle({"stick": {}})
            viewer.zoomTo()

            mol3d_html = viewer.render()
            html_main= viewer._make_html()
            html = f'''
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>3D Protein Visualization</title>
    <script src="https://cdn.tailwindcss.com"></script>
</head>
<body class="font-sans bg-gradient-to-b from-blue-50 to-blue-100 text-blue-800 min-h-screen" style="margin: 0; padding: 0;">
    <!-- Navigation Bar -->
    <nav class="text-white p-4 flex justify-between items-center shadow-lg fixed w-full top-0 z-10" style="background: rgb(37 99 235 / var(--tw-bg-opacity, 1));">
        <a href="#" class="text-3xl font-extrabold">ProteinViz</a>
        <ul class="flex gap-6 text-lg">
            <li><a href="#about" class="hover:underline">About</a></li>
            <li><a href="#services" class="hover:underline">Services</a></li>
            <li><a href="#contact" class="hover:underline">Contact</a></li>
        </ul>
    </nav>

    <!-- Hero Section -->
    <header class="py-24 text-center mt-16">
        <h1 class="text-5xl font-bold text-blue-900">3D Protein Visualization</h1>
        <h2 class="text-2xl mt-4 font-medium text-blue-700">Computational Structural Biology</h2>
        <p class="text-xl mt-4 font-semibold text-blue-800">Visualization of <span class="font-bold">{dump_smile}</span></p>
    </header>

    <!-- 3D View Section -->
    <section id="3dview" class="flex justify-center items-center h-[60vh] p-6">
            {html_main}
       
    </section>

    <!-- Floating Decorative Elements -->
    <div class="absolute top-10 left-5 w-24 h-24 bg-blue-200 rounded-full opacity-50" style="animation: bounce 2s infinite;"></div>
    <div class="absolute bottom-10 right-10 w-20 h-20 bg-blue-200 rounded-full opacity-50" style="animation: pulse 3s infinite;"></div>
    <div class="absolute top-1/3 left-1/4 w-16 h-16 bg-blue-300 rounded-full opacity-40" style="animation: spin 4s linear infinite;"></div>

    <!-- Footer -->
    <footer class="text-white py-6 mt-16 text-center shadow-lg" style="background: rgb(37 99 235 / var(--tw-bg-opacity, 1));">
        <p class="text-lg font-semibold">Developed by <b>N V ANAND SAI KUMAR</b></p>
        <p>&copy; 2023 3D Protein Visualization</p>
        <p><a href="#" class="underline">Privacy Policy</a> | <a href="#" class="underline">Terms of Service</a></p>
    </footer>
</body>
</html>


'''         
            # f = open('templates/main.html',"wt")
            # f.write(html)
            # Redirect to the output page with the visualization
            return html
        else:
            error_message = "Invalid SMILES input. Please try again."
            return render_template('index.html', error_message=error_message)

    return render_template('index.html')

# Define the route for the output page
@app.route('/output')
def output():
    smiles = request.args.get('smiles')
    mol3d_html = request.args.get('mol3d_html')
    return render_template("main.html")

if __name__ == '__main__':
    app.run(debug=True,host='0.0.0.0',port=80)
