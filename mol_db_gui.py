#importing necessary packages
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, Crippen, rdmolops
import sqlite3
import tkinter as tk
from tkinter import ttk, simpledialog
import io
from PIL import Image, ImageTk
import requests
import concurrent.futures
import webbrowser
import csv

# Function to generate and save molecule images
def generate_molecule_image(mol):
    img = Draw.MolToImage(mol)  # Generate image from molecule
    byte_io = io.BytesIO()  # Save image to a byte stream
    img.save(byte_io, format='PNG')  # Save image in PNG format
    byte_io.seek(0)  
    return byte_io.getvalue()  

def get_pubchem_id(mol):
    smiles = Chem.MolToSmiles(mol) # Convert molecule to SMILE
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/TXT" # Obtain PubChem ID using the SMILE
    try:
        response = requests.get(url) # Use PubChem API to get URL
        if response.status_code == 200:
            return response.text.strip()  # Return PubChem CID if found
    except Exception as e:
        return None

# Function to generate images for a batch of molecules in parallel
def generate_molecule_images(mol_data_list):
    with concurrent.futures.ThreadPoolExecutor() as executor:
        images = list(executor.map(generate_molecule_image, [mol_data['mol'] for mol_data in mol_data_list]))
    for idx, image in enumerate(images):
        mol_data_list[idx]['Image'] = image
    return mol_data_list  

# Function to obtain pubchem ids for a batch of molecules in parallel
def fetch_pubchem_ids(mol_data_list):
    pubchem_ids = []
    for mol_data in mol_data_list:
        pubchem_id = get_pubchem_id(mol_data['mol'])
        pubchem_ids.append(pubchem_id)
    return pubchem_ids

'''
def fetch_pubchem_id_parallel(mol_data_list):
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(executor.map(get_pubchem_id, [mol_data['mol'] for mol_data in mol_data_list]))
    return results'''

# Creating function to parse sdf file and extract molecule data from it
def extract_mol_data(mol):
  mol_data = {} # Empty dict
  mol_data['mol'] = mol
  # Extracting and storing molecule properties
  for prop_name in mol.GetPropsAsDict():
    mol_data[prop_name] = mol.GetProp(prop_name)

  # Calculating logP using crippen module of rdkit
  mol_data['LogP'] = Crippen.MolLogP(mol)  #LogP

  # Obtaining hydrogen bond donor and acceptors number using descriptors module of rdkit
  mol_data['HBondDonors'] = Descriptors.NumHDonors(mol)  #Hydrogen Bond Donors
  mol_data['HBondAcceptors'] = Descriptors.NumHAcceptors(mol)  #Hydrogen Bond Acceptors
  mol_data['RotatableBonds'] = Descriptors.NumRotatableBonds(mol) #Number of rotatable bonds
  mol_data['Ring_Count'] = len(mol.GetRingInfo().AtomRings())  # Total number of rings
  mol_data['PSA'] = Descriptors.TPSA(mol)  # Polar Surface Area
  mol_data['FusedRing_Count'] = len(rdmolops.GetSymmSSSR(mol))  # Count of fused rings

  return mol_data #Return the dict

# Creating a function to create the molecule database
def create_mol_db(database):
  create_sql = '''
        --creating a table in database to store the molecule information
        CREATE TABLE IF NOT EXISTS Molecules (
            Image BLOB,
            Name TEXT PRIMARY KEY,
            Mol_ID TEXT,
            Mol_Weight REAL,
            LogP REAL,
            LogD REAL,
            Formula TEXT,
            HBond_Donors INT,
            HBond_Acceptors INT,
            Ring_Count INT,
            PSA REAL,
            FusedRing_Count INT,
            Rotatable_Bonds INT,
            PubChem_ID TEXT            
        );
    ''' # Query to be passed in SQL
  try:
    conn = sqlite3.connect(database) # Connecting to the database
    cur = conn.cursor()
    cur.executescript(create_sql) # Executing above sql command to create database
    conn.commit() # Committing the changes 
    conn.close() # Closing the connection
  except Exception as error: # Incase of an error, informing the user
    print(f'Database could not be created due to : {error}.')
  
def load_mol_data(mol_data_list, database, batch_size=50):
    load_sql = '''
        INSERT INTO Molecules (Image, Name, Mol_ID, Mol_Weight, LogP, LogD, Formula, HBond_Donors, HBond_Acceptors, Ring_Count, PSA, FusedRing_Count, Rotatable_Bonds, PubChem_ID)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    '''
    try:
        conn = sqlite3.connect(database)
        conn.isolation_level = None  # Disable auto-commit for better performance
        cur = conn.cursor()

        # Breaking data into smaller batches for efficiency
        for i in range(0, len(mol_data_list), batch_size):
            batch = mol_data_list[i:i+batch_size]
            batch = generate_molecule_images(batch)
            pubchem_ids = fetch_pubchem_ids(batch)
            data_to_insert = [] # Empty list to store data in batches
            for idx, mol_data in enumerate(batch):
                pubchem_id = pubchem_ids[idx]
                pubchem_link = f"https://pubchem.ncbi.nlm.nih.gov/compound/{pubchem_id}" if pubchem_id else None
                data_to_insert.append((
                    mol_data['Image'], mol_data['Name'].lower(), mol_data['Mol_ID'], round(float(mol_data['MolWeight']), 2),
                    round(float(mol_data['LogP']), 2), round(float(mol_data['LogD']), 2), mol_data['Formula'],
                    mol_data['HBondDonors'], mol_data['HBondAcceptors'], mol_data['Ring_Count'], 
                    round(float(mol_data['PSA']), 2), mol_data['FusedRing_Count'], mol_data['RotatableBonds'], pubchem_link
                ))

            cur.executemany(load_sql, data_to_insert) # Executing multiple queries at once
            conn.commit()  # Commit after each batch

        conn.close()
    except Exception as error:
        print(f'Data could not be loaded into database due to: {error}') # Incase of an error, informing the user
        cur.close()
        conn.close() # Closing the connection

# Creating a fucntion to create the GUI

def create_gui(database):
    app = tk.Tk() 
    app.title("Molecule Database") # Add title
    app.configure(bg="#f8ccf6") #Pick a background colour

    # Style Configuration
    style = ttk.Style()
    style.configure("Treeview.Heading", font=('Georgia', 11, 'bold'), background="#6c0767", foreground="black")
    style.configure("Treeview", font=('Georgia', 10), foreground="black")

    # Create a Treeview widget having custom columns with headings and defined width
    table = ttk.Treeview(app, style="Treeview", columns=('SrNo', 'Name', 'Mol_ID', 'Mol_Weight', 'LogP', 
                                                         'LogD', 'Formula', 'HBond_Donors', 'HBond_Acceptors', 
                                                         'Ring_Count', 'PSA', 'FusedRing_Count', 'Rotatable_Bonds', 'PubChem_ID'))
    table.column('#0', width=0, stretch=tk.NO)
    column_widths = {'SrNo': 50, 'Name': 400, 'Mol_ID': 100, 'Mol_Weight': 100, 'LogP': 100, 
                     'LogD': 100, 'Formula': 120, 'HBond_Donors': 100, 'HBond_Acceptors': 100, 'Ring_Count': 100, 
                     'PSA': 100, 'FusedRing_Count': 100, 'Rotatable_Bonds': 100, 'PubChem_ID': 350}
    for col, width in column_widths.items():
        table.column(col, anchor=tk.W, width=width)

    headings = {'SrNo': 'No.', 'Name': 'Name (Click to view Structure)', 'Mol_ID': 'Mol.ID', 'Mol_Weight': 'Mol.Wt', 
                'LogP': 'LogP', 'LogD': 'LogD', 'Formula': 'Formula', 'HBond_Donors': 'Donors', 
                'HBond_Acceptors': 'Acceptors', 'Ring_Count': 'Rings', 'PSA': 'PSA', 'FusedRing_Count': 'Fused Rings', 
                'Rotatable_Bonds': 'Rot. Bonds', 'PubChem_ID':'PubChem ID (Click to open)'}
    for col, text in headings.items():
        table.heading(col, text=text, anchor=tk.W)

    # Alternating row colors for clarity
    table.tag_configure('oddrow', background='#ffdcf9')
    table.tag_configure('evenrow', background='#f8b0eb')

    # Function to fetch data from database to be loaded into the table
    def fetch_data(query):
        try:
            with sqlite3.connect(database) as conn:
                cur = conn.cursor()
                cur.execute(query)
                rows = cur.fetchall()
                if not rows:
                    print("No data found.")
                return rows
        except Exception as e:
            print("Error fetching data:", e)
            return []
        
    image_refs = {} # Empty dict to store the images

    # Function to open the PubChem URL in a web browser
    def open_pubchem(event, url):
        webbrowser.open(url)

    # Function to pens a new window displaying the molecule image
    def show_molecule_image(mol_id):
        if mol_id in image_refs:
            img_window = tk.Toplevel(app)
            img_window.title("Molecule Structure")
            img_window.geometry("200x200")
            img_label = tk.Label(img_window, image=image_refs[mol_id])
            img_label.pack(expand=True)
    
    # Function to handle clicks on the Treeview to open molecule images or PubChem links
    def on_item_click(event):
        item_id = table.identify_row(event.y)
        col = table.identify_column(event.x)

        if item_id:
            if col == "#2":  # molecule name is in the 2nd column
                show_molecule_image(int(item_id))  # Open image window
            elif col == f"#{len(table['columns'])}":  # Last column for URL
                url = table.item(item_id, "values")[-1]  # Get URL from last column
                open_pubchem(event, url)
    table.bind("<Double-1>", on_item_click)

    displayed_rows = []  # A list to keep track of the currently displayed rows to later export if needed

    # Function to load queries data into the table
    def load_data(query):
        global displayed_rows  
        for row in table.get_children():
            table.delete(row) # Clear the current table
        data = fetch_data(query) # Fetch data based on the query
        displayed_rows = data  # Keep track of the currently displayed data

        # Load data for each row
        for i, row in enumerate(data, 1):
            img_data = row[0]  # BLOB image data
            try:
                img = Image.open(io.BytesIO(img_data))
                img = img.resize((150, 150), Image.LANCZOS)
                img_tk = ImageTk.PhotoImage(img)
                image_refs[i] = img_tk
            except Exception as e:
                print(e)

            table.insert('', 'end', iid=i, values=(i, *row[1:]), 
                        tags=('evenrow' if i % 2 == 0 else 'oddrow',))

    # Function to search molecule by name
    def search_molecule(event=None):
        search_query = search_entry.get().lower()
        # Construct a SQL query with a WHERE clause to filter by the molecule name
        query = f"SELECT * FROM Molecules WHERE Name LIKE '%{search_query}%'"
        load_data(query)

    # Function to show/hide database
    def toggle_table():
        if table.winfo_ismapped():
            table.grid_forget()
        else:
            table.grid(row=2, column=0, columnspan=10, sticky="nsew")
            load_data('SELECT * FROM Molecules')

    # Function to load filtered data
    def filter_data(query):
        load_data(query)

    # Function to obtain user inputed parameters to be used in filters
    def submit_manual_entry():
        mol_weight = float(mol_weight_entry.get()) if mol_weight_entry.get() else 500
        logp = float(logp_entry.get()) if logp_entry.get() else 5
        psa = float(psa_entry.get()) if psa_entry.get() else 200
        hbond_donors = int(hbond_donors_entry.get()) if hbond_donors_entry.get() else 5
        hbond_acceptors = int(hbond_acceptors_entry.get()) if hbond_acceptors_entry.get() else 10
        ring_count = int(ring_count_entry.get()) if ring_count_entry.get() else 4
        fused_ring_count = int(fused_ring_count_entry.get()) if fused_ring_count_entry.get() else 5
        logd_lower = float(logd_lower_entry.get()) if logd_lower_entry.get() else -4
        logd_upper = float(logd_upper_entry.get()) if logd_upper_entry.get() else 4
        rotatable_bonds = int(rotatable_bonds_entry.get()) if rotatable_bonds_entry.get() else 10

        # Modify filter queries dynamically based on the inputs
        global filters
        filters = {
            'Lipinski': f'''SELECT * FROM Molecules WHERE Mol_Weight <= {mol_weight} 
                            AND LogP <= {logp} AND HBond_Donors <= {hbond_donors} 
                            AND HBond_Acceptors <= {hbond_acceptors};''',
            
            'Lead-likeness': f'''SELECT * FROM Molecules WHERE Mol_Weight <= {mol_weight} 
                                AND LogD BETWEEN {logd_lower} AND {logd_upper} 
                                AND Ring_Count <= {ring_count} 
                                AND HBond_Donors <= {hbond_donors} 
                                AND HBond_Acceptors <= {hbond_acceptors} 
                                AND Rotatable_Bonds <= {rotatable_bonds};''',
            
            'Bioavailability': f'''SELECT * FROM Molecules WHERE 
                                (Mol_Weight <= {mol_weight}) + (LogP <= {logp}) + (HBond_Donors <= {hbond_donors}) + 
                                (HBond_Acceptors <= {hbond_acceptors}) + (PSA <= {psa}) + (Rotatable_Bonds <= {rotatable_bonds}) + 
                                (FusedRing_Count <= {fused_ring_count}) >= 6;'''
        }
    
    # Create buttons for the filters
    filter_label = tk.Label(app, text="Filter By:", font=('Georgia', 12, 'bold'), bg="#f8ccf6")
    filter_label.grid(row=1, column=1, padx=10, pady=10)

    filter_button = tk.Button(app, text="Lipinski Rule of 5", command=lambda: filter_data(filters['Lipinski']), bg="#cc1f8f", fg="white", font=('Georgia', 10, 'bold'))
    lead_button = tk.Button(app, text="Lead-likeness", command=lambda: filter_data(filters['Lead-likeness']), bg="#cc1f8f", fg="white", font=('Georgia', 10, 'bold'))
    bio_button = tk.Button(app, text="Bioavailability", command=lambda: filter_data(filters['Bioavailability']), bg="#cc1f8f", fg="white", font=('Georgia', 10, 'bold'))
    filter_button.grid(row=1, column=2, padx=10, pady=10)
    lead_button.grid(row=1, column=3, padx=10, pady=10)
    bio_button.grid(row=1, column=4, padx=10, pady=10)

    # Function to rank data according to different parameters in ascending order
    def rank_data(selected_param):
        query = f'SELECT * FROM Molecules ORDER BY {selected_param} ASC'
        load_data(query)
        
    # Define the ranking options
    rank_options = ['Mol_Weight', 'LogP', 'LogD', 'HBond_Donors', 'HBond_Acceptors', 'Ring_Count', 'PSA', 'FusedRing_Count', 'Rotatable_Bonds']
    
    # Create the dropdown menu for ranking options
    selected_param = tk.StringVar(value="Rank in Ascending order according to")
    rank_dropdown = tk.OptionMenu(app, selected_param, *rank_options, command=rank_data)
    rank_dropdown.config(bg="#cc1f8f", fg="white", font=('Georgia', 10, 'bold'), padx=5, pady=5)

    # Function to export the displayed content of table to a csv or tsv file
    def export_to_file(file_type, filename):
        global displayed_rows
        if displayed_rows:
            with sqlite3.connect(database) as conn:
                cur = conn.cursor()
                cur.execute('PRAGMA table_info(Molecules)') # Get column names from the database schema
                columns = [col[1] for col in cur.fetchall()] # Excuding images
            columns_to_export = columns[1:]
            delimiter = ',' if file_type == 'csv' else '\t' # Determine the delimiter based on the file type
            filename = filename + '.' + file_type # Use the user entered filename

            # Write to the file
            with open(filename, 'w', newline='') as file:
                writer = csv.writer(file, delimiter=delimiter)
                writer.writerow(columns_to_export)
                for row in displayed_rows:
                    writer.writerow(row[1:]) 

    # Function to ask for the filename user input
    def ask_for_filename(file_type):
        filename = simpledialog.askstring("Input", f"Enter the name for the {file_type} file:", parent=app)
        if filename:
            export_to_file(file_type, filename)

    # Function to select csv/tsv from dropdown
    def on_export_choice(selected_format):
        ask_for_filename(selected_format)

    # Export dropdown button
    export_options = ['csv', 'tsv']
    export_selected = tk.StringVar(value="Export Database as")
    export_dropdown = tk.OptionMenu(app, export_selected, *export_options, command=on_export_choice)
    export_dropdown.config(bg="#cc1f8f", fg="white", font=('Georgia', 10, 'bold'), padx=5, pady=5)

    # Function to reset back to original database after filtering or ranking  
    def reset_filter():
        load_data('SELECT * FROM Molecules')
    
    # Create entry boxes for Manual Parameter Entry 
    parameter_label = tk.Label(app, text="Enter Custom Molecule Parameters:", font=('Georgia', 12, 'bold'), bg="#f8ccf6")
    parameter_label.grid(row=3, column=0, padx=10, pady=10, sticky="w")

    # Create entry fields for manual parameter input
    mol_weight_entry = tk.Entry(app, font=('Georgia', 10))
    mol_weight_entry.grid(row=4, column=1, padx=10, pady=5)
    mol_weight_label = tk.Label(app, text="Molecular Weight <=:", font=('Georgia', 10, 'bold'), bg="#f8ccf6")
    mol_weight_label.grid(row=4, column=0, padx=10, pady=5, sticky="w")

    logp_entry = tk.Entry(app, font=('Georgia', 10))
    logp_entry.grid(row=5, column=1, padx=10, pady=5)
    logp_label = tk.Label(app, text="LogP <=:", font=('Georgia', 10, 'bold'), bg="#f8ccf6")
    logp_label.grid(row=5, column=0, padx=10, pady=5, sticky="w")
    
    psa_entry = tk.Entry(app, font=('Georgia', 10))
    psa_entry.grid(row=6, column=1, padx=10, pady=5)
    psa_label = tk.Label(app, text="PSA <=:", font=('Georgia', 10, 'bold'), bg="#f8ccf6")
    psa_label.grid(row=6, column=0, padx=10, pady=5, sticky="w")

    hbond_donors_entry = tk.Entry(app, font=('Georgia', 10))
    hbond_donors_entry.grid(row=7, column=1, padx=10, pady=5)
    hbond_donors_label = tk.Label(app, text="HBond Donors <=:", font=('Georgia', 10, 'bold'), bg="#f8ccf6")
    hbond_donors_label.grid(row=7, column=0, padx=10, pady=5, sticky="w")

    hbond_acceptors_entry = tk.Entry(app, font=('Georgia', 10))
    hbond_acceptors_entry.grid(row=8, column=1, padx=10, pady=5)
    hbond_acceptors_label = tk.Label(app, text="HBond Acceptors <=:", font=('Georgia', 10, 'bold'), bg="#f8ccf6")
    hbond_acceptors_label.grid(row=8, column=0, padx=10, pady=5, sticky="w")

    ring_count_entry = tk.Entry(app, font=('Georgia', 10))
    ring_count_entry.grid(row=4, column=3, padx=10, pady=5)
    ring_count_label = tk.Label(app, text="Ring Count <=:", font=('Georgia', 10, 'bold'), bg="#f8ccf6")
    ring_count_label.grid(row=4, column=2, padx=10, pady=5, sticky="w")

    fused_ring_count_entry = tk.Entry(app, font=('Georgia', 10))
    fused_ring_count_entry.grid(row=5, column=3, padx=10, pady=5)
    fused_ring_count_label = tk.Label(app, text="Fused Ring Count <=:", font=('Georgia', 10, 'bold'), bg="#f8ccf6")
    fused_ring_count_label.grid(row=5, column=2, padx=10, pady=5, sticky="w")

    logd_lower_entry = tk.Entry(app, font=('Georgia', 10))
    logd_lower_entry.grid(row=6, column=3, padx=10, pady=5)
    logd_lower_label = tk.Label(app, text="LogD Lower Value >=:", font=('Georgia', 10, 'bold'), bg="#f8ccf6")
    logd_lower_label.grid(row=6, column=2, padx=10, pady=5, sticky="w")

    logd_upper_entry = tk.Entry(app, font=('Georgia', 10))
    logd_upper_entry.grid(row=7, column=3, padx=10, pady=5)
    logd_upper_label = tk.Label(app, text="LogD Upper Value <=:", font=('Georgia', 10, 'bold'), bg="#f8ccf6")
    logd_upper_label.grid(row=7, column=2, padx=10, pady=5, sticky="w")

    rotatable_bonds_entry = tk.Entry(app, font=('Georgia', 10))
    rotatable_bonds_entry.grid(row=8, column=3, padx=10, pady=5)
    rotatable_bonds_label = tk.Label(app, text="Rotatable Bonds <=:", font=('Georgia', 10, 'bold'), bg="#f8ccf6")
    rotatable_bonds_label.grid(row=8, column=2, padx=10, pady=5, sticky="w")

    # Button to submit the entries
    submit_button = tk.Button(app, text="Submit Custom Parameters", command=submit_manual_entry, bg="#cc1f8f", fg="white", font=('Georgia', 10, 'bold'))
    submit_button.grid(row=9, column=1, padx=10, pady=10)

    # Button to search by molecule name 
    search_label = tk.Label(app, text="Search by Molecule Name:", font=('Georgia', 12, 'bold'), bg="#f8ccf6")
    search_entry = tk.Entry(app, font=('Georgia', 10))
    search_entry.bind("<KeyRelease>", search_molecule)  # Trigger search when the user types

    # Buttons for showing/hiding database and reseting
    show_button = tk.Button(app, text="Show/Hide Database", command=toggle_table, bg="#cc1f8f", fg="white", font=('Georgia', 10, 'bold'))
    reset_button = tk.Button(app, text="Reset", command=reset_filter, bg="#cc1f8f", fg="white", font=('Georgia', 10, 'bold'))
    
    #Arranging remaining buttons in first 2 rows
    show_button.grid(row=0, column=1, padx=10, pady=10)
    reset_button.grid(row=0, column=2, padx=10, pady=10)
    export_dropdown.grid(row=0, column=3, padx=10, pady=10)
    rank_dropdown.grid(row=0, column=4, padx=10, pady=10)
    search_label.grid(row=0, column=0, padx=10, pady=10, sticky="w")
    search_entry.grid(row=1, column=0, padx=10, pady=10, sticky="w")

    # Adjust the grid configurations for resizing
    app.grid_rowconfigure(2, weight=1)  # Allow the table area to expand
    app.grid_columnconfigure(0, weight=1)  # Allow the buttons area to expand
    app.grid_columnconfigure(5, weight=2)  # Add weight for rank dropdown to stretch properly

    # Scrollbars on both axes
    vsb = ttk.Scrollbar(app, orient="vertical", command=table.yview)
    vsb.grid(row=2, column=11, sticky="ns")
    table.configure(yscrollcommand=vsb.set)

    hsb = ttk.Scrollbar(app, orient="horizontal", command=table.xview)
    hsb.grid(row=3, column=0, columnspan=10, sticky="ew")
    table.configure(xscrollcommand=hsb.set)

    app.mainloop()


# Loading the SDF file and extracting data from it to be stored into the database
'''database = 'molecules20.db'

create_mol_db(database)

with Chem.SDMolSupplier('Molecules20.sdf') as mols:
  mol_data_list = []
  for mol in mols:
        if mol is not None:
          mol_data = extract_mol_data(mol)
          mol_data_list.append(mol_data)
          
load_mol_data(mol_data_list, database)'''

# Displaying the created database in GUI
create_gui('molecules20.db')

