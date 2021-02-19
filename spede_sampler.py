import PySimpleGUI as sg
import os
import random
import shutil
from Bio import AlignIO
from Bio.Phylo.Applications import RaxmlCommandline
from Bio.Phylo.Applications import FastTreeCommandline

sg.theme('DarkBlue 4')

global ingroup_display
ingroup_display = ""
global outgroup_display
outgroup_display = ""
global iterations
iterations = 0
global no_seqs
no_seqs = 0
global content
content = []
global content2
content2 = []

######################################################################################################################
# MENU OPTIONS
######################################################################################################################

menu_def = [
           ['Resample Fasta files', ['Open Ingroup Fasta File', 'Print Ingroup Fasta File to Screen',
            'Open Outgroup Fasta File', 'Print Outgroup Fasta File to Screen', 'Nucleotide Composition',
                                     'Run Resampling', 'Fasta Resampling Help']],
           ['RAxML', ['Convert Fasta to Phylip', 'Run RAxML', 'RAxML Help']],
           ['FastTree', ['Run FastTree', 'FastTree Help']],
           ['GMYC', ['Run GMYC']],
           ['About', ['Overview', 'Citation', 'Contact']], ]

######################################################################################################################
# MAIN WINDOW LAYOUT
######################################################################################################################

layout_main = [
    [sg.Menu(menu_def, tearoff=False)],
    [sg.Image('spede_sampler_image.png')]
]

######################################################################################################################
# CREATE MAIN WINDOW
######################################################################################################################

window_main = sg.Window('SPEDE-SAMPLER', layout_main, size=(560, 350), icon='bug.ico')

while True:
    event, values = window_main.read()

    if event is None:
        break

######################################################################################################################
# HELP WITH FASTA FILE INPUT FOR RESAMPLING
######################################################################################################################

    if event == 'Overview':
        overview = open('overview.txt', 'r')
        overview_content = overview.read()
        overview.close()
        sg.popup_scrolled(overview_content, title='Fasta Resampling Help', size=(100, 25))

######################################################################################################################
# HELP WITH FASTA FILE INPUT FOR RESAMPLING
######################################################################################################################

    if event == 'Fasta Resampling Help':
        help_fasta_resampling = open('help_fasta_resampling.txt', 'r')
        help_fasta_resampling_content = help_fasta_resampling.read()
        help_fasta_resampling.close()
        sg.popup_scrolled(help_fasta_resampling_content, title='Fasta Resampling Help', size=(100, 25))

######################################################################################################################
# HELP WITH RAxML
######################################################################################################################

    if event == 'RAxML Help':
        help_raxml = open('help_raxml.txt', 'r')
        help_raxml_content = help_raxml.read()
        help_raxml.close()
        sg.popup_scrolled(help_raxml_content, title='RAxML Help', size=(100, 25))

######################################################################################################################
# HELP WITH FASTTREE
######################################################################################################################

    if event == 'FastTree Help':
        help_fasttree = open('help_fasttree.txt', 'r')
        help_fasttree_content = help_fasttree.read()
        help_fasttree.close()
        sg.popup_scrolled(help_fasttree_content, title='FastTree Help', size=(100, 25))

######################################################################################################################
# HELP WITH FASTTREE
######################################################################################################################

    if event == 'Run GMYC':
        help_gmyc = open('help_gmyc.txt', 'r')
        help_gmyc_content = help_gmyc.read()
        help_gmyc.close()
        sg.popup_scrolled(help_gmyc_content, title='GMYC Help', size=(100, 25))

######################################################################################################################
# READ IN AN INGROUP FASTA FILE
######################################################################################################################

    if event == 'Open Ingroup Fasta File':
        # use this popup option for menu items, not for a button like FolderBrowse
        ingroup_fasta_filename = sg.popup_get_file('file to open', no_window=True, file_types=((".fasta", ".fasta"),
                                                                                               (".fas", ".fas")))
        if len(ingroup_fasta_filename) > 0:
            ingroup_fasta_file = open(ingroup_fasta_filename, 'r+')
            sg.popup_notify('Ingroup Fasta file ' + ingroup_fasta_filename + ' successfully read in',
                            location=(600, 400), display_duration_in_ms=4000,
                            fade_in_duration=10, alpha=1)
            ingroup_display = ingroup_fasta_file.read()
            content = ingroup_display.splitlines()
            ingroup_fasta_file.close()
            content2 = []

######################################################################################################################
# PRINT INGROUP FASTA FILE TO THE SCREEN
######################################################################################################################

    if event == 'Print Ingroup Fasta File to Screen':
        if len(ingroup_display) == 0:
            sg.popup_ok('File Upload Required', 'Please upload an ingroup Fasta file to display', icon='bug.ico')
        else:
            sg.popup_scrolled(ingroup_display, title=ingroup_fasta_filename, size=(100, 100))

######################################################################################################################
# READ IN AN OUTGROUP FASTA FILE
######################################################################################################################

    if event == 'Open Outgroup Fasta File':
        outgroup_fasta_filename = sg.popup_get_file('file to open', no_window=True, file_types=((".fasta", ".fasta"),
                                                                                                (".fas", ".fas")))
        if len(outgroup_fasta_filename) > 0:
            outgroup_fasta_file = open(outgroup_fasta_filename, 'r+')
            sg.popup_notify('Outgroup Fasta file ' + outgroup_fasta_filename + ' successfully read in',
                            location=(600, 400), display_duration_in_ms=4000,
                            fade_in_duration=10, alpha=1)
            outgroup_display = outgroup_fasta_file.read()

######################################################################################################################
# PRINT OUTGROUP FASTA FILE TO THE SCREEN
######################################################################################################################

    if event == 'Print Outgroup Fasta File to Screen':
        if len(outgroup_display) == 0:
            sg.popup_ok('File Upload Required', 'Please upload an outgroup Fasta file to display', icon='bug.ico')
        else:
            # noinspection PyUnboundLocalVariable
            sg.popup_scrolled(outgroup_display, title=outgroup_fasta_filename, size=(100, 100))

######################################################################################################################
# NUCLEOTIDE COMPOSITION
######################################################################################################################

    if event == 'Nucleotide Composition':

        if len(ingroup_display) == 0:
            sg.popup_ok('File Upload Required', 'Please upload an ingroup Fasta file before proceeding.',
                        icon='bug.ico')
        else:
            ###########################################
            # Count the number of sequences in the file
            ###########################################
            n = 0
            for line in content:
                if line.startswith(">"):
                    n += 1

            ############################################
            # Count the base numbers for each nucleotide
            ############################################

            j = 0
            seq_length = []
            seq_name = []

            adenosine = 0
            cytosine = 0
            guanine = 0
            thymine = 0
            ambiguous = 0
            gap = 0

            # This loop keeps a count of each character in the sequence (A, C, T, G, - and N)
            for line in content:
                j += 1
                # this makes sure that only even row lengths are stored (omitting the line with the sequence name):
                # store each sequence name
                if j % 2 != 0:
                    seq_name.append(line)
                if j % 2 == 0:
                    seq_length.append(len(line))  # append the length of the line to the object "seq_length"
                    # Now go through each character in the sequence, and check its base
                    for i in range(0, len(line)):
                        if line[i] == "A" or line[i] == "a":
                            adenosine += 1
                        elif line[i] == "C" or line[i] == "c":
                            cytosine += 1
                        elif line[i] == "G" or line[i] == "g":
                            guanine += 1
                        elif line[i] == "T" or line[i] == "t":
                            thymine += 1
                        elif line[i] == "N" or line[i] == "n" or line[i] == "?" or line[i] == "R" or line[i] == "Y" or \
                                line[i] == "S" or line[i] == "W" or line[i] == "K" or line[i] == "M" or \
                                line[i] == "B" or line[i] == "D" or line[i] == "H" or line[i] == "V" or \
                                line[i] == "r" or line[i] == "y" or line[i] == "s" or line[i] == "w" or \
                                line[i] == "k" or line[i] == "m" or line[i] == "b" or line[i] == "d" or \
                                line[i] == "h" or line[i] == "v":
                            ambiguous += 1
                        elif line[i] == "-":
                            gap += 1

            sum_nucs = adenosine + cytosine + guanine + thymine
            percent_ambig = ambiguous / (n * seq_length[0]) * 100
            percent_gap = gap / (n * seq_length[0]) * 100

            ############################################
            # Check if sequence lengths are equal
            ############################################

            first_ele = seq_length[0]  # compare all sequences to the first one
            chk = True

            for item in seq_length:
                if first_ele != item:
                    chk = False
                    break
            if chk:
                length_check = 'All sequences are equal in length, with ' + str(seq_length[0]) + ' nucleotide bases.'
            else:
                # join the two lists containing sequence names and their lengths so that they can display side by side
                output_lens = "\n".join("{} {}".format(x, y) for x, y in zip(seq_name, seq_length))
                length_check = '*NOTE* Sequences are not equal in length. Please ensure that they are, ' \
                               'and re-upload the file. See below for the lengths of each sequence.' + \
                               '\n' + '\n' + output_lens

            ############################################
            # Write the results to the screen
            ############################################

            sg.popup_scrolled('There are ' + str(n) + ' sequences in the ingroup FASTA file. \n\nNucleotide base '
                                                      'compositions for the ingroup sequences are: \n\n(A) ' +
                              str(round(adenosine / sum_nucs * 100)) + '%' + '\n(C) ' +
                              str(round(cytosine / sum_nucs * 100)) + '%' + '\n(G) ' +
                              str(round(guanine / sum_nucs * 100)) + '%' + '\n(T) ' +
                              str(round(thymine / sum_nucs * 100)) + '%' + '\nAmbiguous (N) ' +
                              str(round(percent_ambig)) + '%,' + '\nGaps (- or ?) ' + str(round(percent_gap)) +
                              '%' + '\n', length_check, title='Summary Information', size=(100, 100))

######################################################################################################################
# FASTA SEQUENCE RESAMPLING
######################################################################################################################

    if event == 'Run Resampling':

        if len(ingroup_display) == 0:
            sg.popup_ok('File Upload Required', 'Please upload an ingroup Fasta file before running the analysis.',
                        icon='bug.ico')

######################################################################################################################
# WINDOW LAYOUT FOR FASTA RESAMPLING
######################################################################################################################
        
        else:
            layout_fasta_resample = [
                [sg.Text('Number of sequences to resample:')],
                [sg.InputText(key='seqs', size=(10, 10))],
                [sg.Text('Number of iterations:')],
                [sg.InputText(key='iter', size=(10, 10))],
                [sg.Checkbox('Append Outgroup Sequences?', key='inc_outgroups', default=False)],
                [sg.Input(), sg.FolderBrowse('Browse')],
                [sg.Submit('Run')],
                [sg.ProgressBar(iterations, orientation='h', size=(20, 20), key='progbar_fasta')]
            ]
            # open a new window with settings
            window_fasta_resample = sg.Window('Resample Fasta Sequences', layout_fasta_resample, size=(500, 250))
            while True:
                event, values = window_fasta_resample.read()
                if event is None:
                    break

                # get the input values from the user
                no_seqs = values['seqs']
                iterations = values['iter']

######################################################################################################################
# RUN THE RESAMPLING ANALYSIS
######################################################################################################################

                if event == 'Run':

                    if values['inc_outgroups'] and len(outgroup_display) == 0:
                        sg.popup_ok('Outgroups?', 'You have indicated that outgroups should be included, '
                                                  'but you have not uploaded a Fasta file for them.')

                    elif len(no_seqs) == 0 or len(iterations) == 0:
                        sg.popup_ok('Iterations?', 'You have not specified values for the number of '
                                                   'sequences and iterations.')

                    elif int(no_seqs) <= 0 or int(iterations) <= 0:
                        sg.popup_ok('Negative Values', 'You have inserted negative numbers.')

                    else:
                        save_dir = values['Browse']
                        os.chdir(save_dir)
                        # This if statement checks if the Iterations folder already exists, and then overwrites it
                        # if it does
                        if os.path.exists(save_dir + '\\' + 'Iterations/'):
                            shutil.rmtree(save_dir + '\\' + 'Iterations/')
                        os.mkdir('Iterations/')

                        for i in range(1, len(content) + 1):
                            if i % 2 == 0:
                                content2.append(content[i - 2].strip() + " " + content[i - 1].strip())

                        f1 = ()
                        counter = 0

                        while counter < int(iterations):

                            r_shuffle = random.sample(content2, k=int(no_seqs))

                            for line in r_shuffle:
                                # splitting the line here takes the fasta file back to its original form,
                                # with the name on one line, and the sequence below it
                                y = line.split()
                                f1 = open("Iterations/iteration_" + str(counter + 1) + ".fas", "a")
                                for element in y:
                                    f1.write(element + "\n")
                            # if outgroups is selected, append them to the end of each file
                            if values['inc_outgroups']:
                                outgroups = outgroup_display.splitlines()
                                for lines in outgroups:
                                    f1.write(lines + "\n")
                                # noinspection PyUnboundLocalVariable
                                outgroup_fasta_file.close()
                            counter += 1
                            f1.close()

                            window_fasta_resample['progbar_fasta'].update_bar(counter + 1)

                        sg.popup_notify(str(iterations) + ' files successfully written to ' + save_dir +
                                        '/Iterations', location=(650, 400), display_duration_in_ms=6000,
                                        fade_in_duration=10, alpha=1)

######################################################################################################################
# CONVERT FASTA TO PHYLIP FOR RAXML
######################################################################################################################

    if event == 'Convert Fasta to Phylip':
        convert_folder = sg.popup_get_folder('fasta files to convert', no_window=True)

        if len(convert_folder) > 0:
            full_filenames = []
            short_filename = []

            for file in os.listdir(convert_folder):
                name = os.fsdecode(file)
                if name.endswith((".fas", ".fasta", ".FASTA")):
                    full_filenames.append(convert_folder + "\\" + name)
                    short_filename.append(name)

            if len(full_filenames) == 0:
                sg.popup_ok('No .fas files', "There are no Fasta files in the selected folder.")

######################################################################################################################
# start of try-catch block to avoid the input of sequences of unequal length or repeated names
######################################################################################################################
            try:
                os.chdir(convert_folder)
                if os.path.exists(convert_folder + '\\' + 'Phylip_conversions/'):
                    shutil.rmtree(convert_folder + '\\' + 'Phylip_conversions/')
                os.mkdir('Phylip_conversions/')

                i = 0
                for f in full_filenames:
                    if f.endswith(".fas"):
                        AlignIO.convert(f, "fasta", "Phylip_conversions/" + short_filename[i].replace(".fas", "")
                                        + ".phy", "phylip", "DNA")
                    elif f.endswith(".fasta"):
                        AlignIO.convert(f, "fasta",
                                        "Phylip_conversions/" + short_filename[i].replace(".fasta", "") + ".phy",
                                        "phylip", "DNA")
                    elif f.endswith(".FASTA"):
                        AlignIO.convert(f, "fasta",
                                        "Phylip_conversions/" + short_filename[i].replace(".FASTA", "") + ".phy",
                                        "phylip", "DNA")
                    i += 1

                sg.popup_notify(str(len(full_filenames)) + " successfully converted from .fas to .phy files. "
                                                           "See the \'Phylip_conversions\' folder in the " +
                                convert_folder + ' directory.', location=(650, 400), display_duration_in_ms=8000,
                                fade_in_duration=10, alpha=1)

            except ValueError:
                # sequence names must not be more than 10 chars, otherwise the rest is cut off, and then there may be
                # repeated names
                sg.popup_ok('File Issue', 'Unequal sequence lengths or repeated names. '
                                          'Please check your file and re-upload.')

######################################################################################################################
# end of try-catch block
######################################################################################################################

######################################################################################################################
# RUN RAxML
# user must make sure that the folder with phylip files is in the same directory as this program in order
# for the raxmlHPC.exe file to be accessed
######################################################################################################################

    if event == 'Run RAxML':
        layout_raxml = [
            [sg.Input(), sg.FolderBrowse(key='raxml_infolder')],
            [sg.Text('Model:')],
            [sg.InputCombo(('GTRCAT', 'GTRGAMMA'),
                           size=(20, 1), key='raxml_model')],
            [sg.Text('Number of bootstrap iterations:')],
            [sg.InputText(key='bootstraps_raxml', size=(10, 10))],
            [sg.Submit('Run', key='run_raxml')],
            [sg.ProgressBar(10, orientation='h', size=(20, 20), key='progbar_raxml')]
        ]
        # open a new window with settings
        window_raxml = sg.Window('RAxML', layout_raxml, size=(500, 250))

        while True:
            event, values = window_raxml.read()
            if event is None:
                break

            raxml_input_folder = values['raxml_infolder']
            boots_raxml = values['bootstraps_raxml']
            model = values['raxml_model']

            if event == 'run_raxml':

                try:
                    boots_raxml_int = int(boots_raxml)
                except ValueError:
                    sg.popup_ok('Error', 'Please input an integer for bootstrap number.')
                    continue

                if len(raxml_input_folder) == 0:
                    sg.popup_ok('Select Folder', 'Please select a folder containing your Phylip files.')
                elif len(model) == 0:
                    sg.popup_ok('Model Selection', 'Please select an evolutionary model from the dropdown menu.')
                elif len(boots_raxml) == 0:
                    sg.popup_ok('Bootstrap Value', 'Number of bootstraps missing. Please check.')
                elif int(boots_raxml) < 0:
                    sg.popup_ok('Negative Value', 'Please insert a positive bootstrap value.')

                else:
                    full_filenames = []
                    short_filenames = []

                    if os.path.exists(raxml_input_folder + '\\' + 'RAxML_output/'):
                        shutil.rmtree(raxml_input_folder + '\\' + 'RAxML_output/')

                    os.mkdir(raxml_input_folder + '\\' + 'RAxML_output/')

                    for file in os.listdir(raxml_input_folder):
                        name = os.fsdecode(file)
                        if name.endswith(".phy"):
                            full_filenames.append(raxml_input_folder + "\\" + name)
                            short_filenames.append(name)

                    if len(full_filenames) == 0:
                        sg.popup_ok("There are no .phy files in the selected folder.")

                    else:
                        try:
                            j = 0
                            # Loop through each file in the desired folder
                            for file in full_filenames:
                                raxml_input = file
                                # specify the name of the output tree file:
                                raxml_outfile = short_filenames[j].replace(".phy", "")
                                cmd_raxml = RaxmlCommandline(sequences=raxml_input, name=raxml_outfile, model=model,
                                                             num_replicates=int(boots_raxml),
                                                             working_dir=raxml_input_folder + '\\' + 'RAxML_output/')
                                cmd_raxml()
                                j += 1
                                window_raxml['progbar_raxml'].update_bar(j + 1)

                            sg.popup_ok('Success!', 'ML trees successfully created using ' + str(boots_raxml) +
                                        ' bootstrap repeats. These have been written to the ' + raxml_input_folder +
                                        '\\' + 'RAxML_output folder.')

                        except:
                            sg.popup_ok('RAxML encountered an error. This is most likely due to your input '
                                        'files having too few sequences, or very short sequence lengths.')

######################################################################################################################
# FASTTREE
######################################################################################################################

    if event == 'Run FastTree':
        layout_fasttree = [
            [sg.Input(), sg.FolderBrowse(key='fasttree_infolder')],
            [sg.Text('Number of bootstrap iterations:')],
            [sg.InputText(key='bootstraps_fasttree', size=(10, 10))],
            [sg.Checkbox('GTR', key='gtr'), sg.Checkbox('Gamma', key='gamma'),
             sg.Checkbox('Nucleotide', default=True, key='nt', tooltip='Deselect for protein data')],
            [sg.Submit('Run', key='run_fasttree')],
            [sg.ProgressBar(10, orientation='h', size=(20, 20), key='progbar_fasttree')]
        ]
        # open a new window with settings
        window_fasttree = sg.Window('FastTree', layout_fasttree, size=(500, 250))

        while True:
            event, values = window_fasttree.read()
            if event is None:
                break

            fasttree_input_folder = values['fasttree_infolder']
            boots_fasttree = values['bootstraps_fasttree']

            if event == 'run_fasttree':

                try:
                    boots_fasttree_int = int(boots_fasttree)
                except ValueError:
                    sg.popup_ok('Error', 'Bootstrap value must be an integer.')
                    continue

                if len(fasttree_input_folder) == 0:
                    sg.popup_ok('Select Folder', 'Please select a folder containing your Fasta files.')
                elif len(boots_fasttree) == 0:
                    sg.popup_ok('Bootstrap Value', 'Please input the number of bootstraps you wish to run.')
                elif int(boots_fasttree) < 0:
                    sg.popup_ok('Negative Value', 'Please insert a positive bootstrap value.')
                else:
                    full_filenames = []
                    short_filenames = []

                    for file in os.listdir(fasttree_input_folder):
                        name = os.fsdecode(file)
                        if name.endswith(".fas") or name.endswith(".fasta"):
                            full_filenames.append(fasttree_input_folder + "\\" + name)
                            short_filenames.append(name)

                    if len(full_filenames) == 0:
                        sg.popup_ok("There are no Fasta files in the selected folder.")

                    else:
                        # create an output folder to save trees to
                        if os.path.exists(fasttree_input_folder + '\\' + 'FastTree_output/'):
                            shutil.rmtree(fasttree_input_folder + '\\' + 'FastTree_output/')

                        os.mkdir(fasttree_input_folder + '\\' + 'FastTree_output/')

                        # set the file path for the FastTree.exe file
                        fasttree_exe = os.path.dirname(os.path.abspath('spede_sampler.py')) + '\\' + 'FastTree.exe'

                    try:
                        q = 0
                        # Loop through each file in the desired folder
                        for file in full_filenames:
                            fasttree_fastafile = file
                            # specify the location of the output tree file:
                            if file.endswith('.fas'):
                                fasttree_outfile = fasttree_input_folder + '\\' + 'FastTree_output/' + \
                                                    short_filenames[q].replace('.fas', '.tre')
                            elif file.endswith('.FASTA'):
                                fasttree_outfile = fasttree_input_folder + '\\' + 'FastTree_output/' + \
                                                    short_filenames[q].replace('.FASTA', '.tre')
                            elif file.endswith('.fasta'):
                                fasttree_outfile = fasttree_input_folder + '\\' + 'FastTree_output/' + \
                                                    short_filenames[q].replace('.fasta', '.tre')
                            # Modify parameters as desired
                            cmd_fasttree = FastTreeCommandline(fasttree_exe, nt=values['nt'], gtr=values['gtr'],
                                                               gamma=values['gamma'], input=fasttree_fastafile,
                                                               out=fasttree_outfile, boot=int(boots_fasttree))
                            cmd_fasttree()
                            q += 1
                            window_fasttree['progbar_fasttree'].update_bar(q + 1)

                        sg.popup_ok('Success!', 'ML trees successfully created with ' + str(boots_fasttree) +
                                    ' bootstrap repeats. These have been written to the ' + fasttree_input_folder +
                                    '\\' + 'FastTree_output folder.')

                    except:
                        sg.popup_ok('FastTree encountered an error.')

######################################################################################################################
# CONTACT DETAILS WINDOW
######################################################################################################################

    if event == 'Contact':
        sg.popup_ok('Report bugs or make suggestions',  'Clarke van Steenderen \nemail: vsteenderen@gmail.com or '
                                                        'g14v1511@campus.ru.ac.za')

######################################################################################################################
# CITATION DETAILS WINDOW
######################################################################################################################

    if event == 'Citation':
        sg.popup_ok('Cite', 'Please cite this program as follows: \nvan Steenderen, CJM. SPEDE-SAMPLER: '
                            'assess sampling effects on species delimitation, version 1.1, 2020.')

window_main.close()
