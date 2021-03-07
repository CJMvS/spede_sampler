import PySimpleGUI as sg
import os
import random
import shutil
from Bio import AlignIO
from Bio.Phylo.Applications import RaxmlCommandline
from Bio.Phylo.Applications import FastTreeCommandline

sg.theme('Reddit')

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
global n  # number of sequences in the uploaded fasta file
global short_filenames_raxml
short_filenames_raxml = []
global short_filenames_fasttree
short_filenames_fasttree =[]
global full_filenames_raxml
full_filenames_raxml = []
global full_filenames_fasttree
full_filenames_fasttree = []

# set the file path for the FastTree.exe file
global fasttree_exe
fasttree_exe = os.path.dirname(os.path.abspath('spede_sampler.py')) + '\\' + 'FastTree.exe'
global raxmlHPC_path
raxmlHPC_path = os.path.dirname(os.path.abspath('spede_sampler.py'))

######################################################################################################################
# MENU OPTIONS
######################################################################################################################

menu_def = [
           ['Resample Fasta files', ['Open Ingroup Fasta File', 'Print Ingroup Fasta File to Screen',
            'Open Outgroup Fasta File', 'Print Outgroup Fasta File to Screen', 'Nucleotide Composition',
                                     'Run Resampling', 'Fasta Resampling Help']],
           ['RAxML', ['Convert Fasta to Phylip', 'Run RAxML (single Phylip file)', 'Run RAxML (multiple Phylip files)',
                      'RAxML Help']],
           ['FastTree', ['Run FastTree (single Fasta file)', 'Run FastTree (multiple Fasta files)', 'FastTree Help']],
           ['GMYC', ['Run GMYC']],
           ['About', ['Overview', 'Citation', 'Contact']], ]

######################################################################################################################
# MAIN WINDOW LAYOUT
######################################################################################################################

layout_main = [
    [sg.Menu(menu_def, tearoff=False, background_color="lightgreen")],
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
# OVERVIEW
######################################################################################################################

    if event == 'Overview':
        with open('overview.txt', 'r') as overview:
            overview_content = overview.read()
        sg.popup_scrolled(overview_content, title='Overview of SPEDE-SAMPLER', size=(100, 25))

######################################################################################################################
# HELP WITH FASTA FILE INPUT FOR RESAMPLING
######################################################################################################################

    if event == 'Fasta Resampling Help':
        with open('help_fasta_resampling.txt', 'r') as help_fasta_resampling:
            help_fasta_resampling_content = help_fasta_resampling.read()
        sg.popup_scrolled(help_fasta_resampling_content, title='Fasta Resampling Help', size=(100, 25))

######################################################################################################################
# HELP WITH RAxML
######################################################################################################################

    if event == 'RAxML Help':
        with open('help_raxml.txt', 'r') as help_raxml:
            help_raxml_content = help_raxml.read()
        sg.popup_scrolled(help_raxml_content, title='RAxML Help', size=(100, 25))

######################################################################################################################
# HELP WITH FASTTREE
######################################################################################################################

    if event == 'FastTree Help':
        with open('help_fasttree.txt', 'r') as help_fasttree:
            help_fasttree_content = help_fasttree.read()
        sg.popup_scrolled(help_fasttree_content, title='FastTree Help', size=(100, 25))

######################################################################################################################
# HELP WITH GMYC
######################################################################################################################

    if event == 'Run GMYC':
        with open('help_gmyc.txt', 'r') as help_gmyc:
            help_gmyc_content = help_gmyc.read()
        sg.popup_scrolled(help_gmyc_content, title='GMYC Help', size=(100, 25))

######################################################################################################################
# READ IN AN INGROUP FASTA FILE
######################################################################################################################

    if event == 'Open Ingroup Fasta File':
        # use this popup option for menu items, not for a button like FolderBrowse
        ingroup_fasta_filename = sg.popup_get_file('file to open', no_window=True, file_types=((".fasta", ".fasta"),
                                                                                               (".fas", ".fas")))
        if len(ingroup_fasta_filename) > 0:

            with open(ingroup_fasta_filename, 'r+') as ingroup_fasta_file:
                ingroup_display = ingroup_fasta_file.read()

            content = ingroup_display.splitlines()
            content2 = []

            for i in range(1, len(content) + 1):
                if i % 2 == 0:
                    content2.append(content[i - 2].strip() + " " + content[i - 1].strip())

            # get the number of sequences in the file
            n = 0
            for line in content:
                if line.startswith(">"):
                    n += 1

            if n == 0:
                warning_fasta_len = '*Warning* Your file is empty. Please check this before continuing.'
            else:
                warning_fasta_len = ''

            sg.popup_ok('File successfully uploaded!',
                        'Ingroup Fasta file {0} successfully read in. The file contains {1} sequences.\n\n{2}'
                        .format(os.path.basename(ingroup_fasta_filename), str(n), warning_fasta_len))

######################################################################################################################
# PRINT INGROUP FASTA FILE TO THE SCREEN
######################################################################################################################

    if event == 'Print Ingroup Fasta File to Screen':
        if len(ingroup_display) == 0:
            sg.popup_ok('File Issue', 'You have either not uploaded an ingroup Fasta file to display, or the Fasta '
                                      'file is empty.', icon='bug.ico')
        else:
            sg.popup_scrolled(ingroup_display, title=os.path.basename(ingroup_fasta_filename), size=(100, 150))

######################################################################################################################
# READ IN AN OUTGROUP FASTA FILE
######################################################################################################################

    if event == 'Open Outgroup Fasta File':
        outgroup_fasta_filename = sg.popup_get_file('file to open', no_window=True, file_types=((".fasta", ".fasta"),
                                                                                                (".fas", ".fas")))
        if len(outgroup_fasta_filename) > 0:
            with open(outgroup_fasta_filename, 'r+') as outgroup_fasta_file:
                outgroup_display = outgroup_fasta_file.read()
            sg.popup_ok('Success!', 'Outgroup Fasta file ' + os.path.basename(outgroup_fasta_filename) +
                        ' successfully read in')

######################################################################################################################
# PRINT OUTGROUP FASTA FILE TO THE SCREEN
######################################################################################################################

    if event == 'Print Outgroup Fasta File to Screen':
        if len(outgroup_display) == 0:
            sg.popup_ok('File Upload Required', 'Please upload an outgroup Fasta file to display', icon='bug.ico')
        else:
            # noinspection PyUnboundLocalVariable
            sg.popup_scrolled(outgroup_display, title=os.path.basename(outgroup_fasta_filename), size=(100, 100))

######################################################################################################################
# NUCLEOTIDE COMPOSITION
######################################################################################################################

    if event == 'Nucleotide Composition':

        if len(ingroup_display) == 0:
            sg.popup_ok('File Issue', 'You have either not uploaded an ingroup Fasta file to display, or the Fasta '
                                      'file is empty.', icon='bug.ico')
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
                        if line[i] == "C" or line[i] == "c":
                            cytosine += 1
                        if line[i] == "G" or line[i] == "g":
                            guanine += 1
                        if line[i] == "T" or line[i] == "t":
                            thymine += 1
                        if line[i] == "N" or line[i] == "n" or line[i] == "?" or line[i] == "R" or line[i] == "Y" or \
                           line[i] == "S" or line[i] == "W" or line[i] == "K" or line[i] == "M" or \
                           line[i] == "B" or line[i] == "D" or line[i] == "H" or line[i] == "V" or \
                           line[i] == "r" or line[i] == "y" or line[i] == "s" or line[i] == "w" or \
                           line[i] == "k" or line[i] == "m" or line[i] == "b" or line[i] == "d" or \
                           line[i] == "h" or line[i] == "v":
                            ambiguous += 1
                        if line[i] == "-":
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
                              str(round(adenosine / sum_nucs * 100, 2)) + '%' + '\n(C) ' +
                              str(round(cytosine / sum_nucs * 100, 2)) + '%' + '\n(G) ' +
                              str(round(guanine / sum_nucs * 100, 2)) + '%' + '\n(T) ' +
                              str(round(thymine / sum_nucs * 100, 2)) + '%' + '\nAmbiguous (N) ' +
                              str(round(percent_ambig, 2)) + '%,' + '\nGaps (- or ?) ' + str(round(percent_gap, 2)) +
                              '%' + '\n', length_check, title='Summary Information', size=(100, 100))

######################################################################################################################
# FASTA SEQUENCE RESAMPLING
######################################################################################################################

    if event == 'Run Resampling':

        if len(ingroup_display) == 0:
            sg.popup_ok('File Issue', 'You have either not uploaded an ingroup Fasta file to display, or the Fasta '
                                      'file is empty.', icon='bug.ico')

######################################################################################################################
# WINDOW LAYOUT FOR FASTA RESAMPLING
######################################################################################################################

        else:
            layout_fasta_resample = [
                [sg.Text('Percentage of the data to resample of the ' + str(n) + ' sequences inputted:')],
                [sg.InputText(key='seqs', size=(10, 10)), sg.Text('%')],
                [sg.Text('Number of iterations:')],
                [sg.InputText(key='iter', size=(10, 10))],
                [sg.Checkbox('Append Outgroup Sequences?', key='inc_outgroups', default=False)],
                [sg.Checkbox('Set a random seed?', key='set_seed', default=False, tooltip='Select this if you want '
                            'the resampling results to be reproducible if you run this again.')],
                [sg.Text('Select a Folder to which Results are Written:')],
                [sg.Input(), sg.FolderBrowse('Browse')],
                [sg.Submit('Run')],
                [sg.ProgressBar(n, orientation='h', size=(20, 20), key='progbar_fasta')]
            ]
            # open a new window with settings
            window_fasta_resample = sg.Window('Resample Fasta Sequences', layout_fasta_resample, size=(500, 300))
            while True:
                event_resample, values_resample = window_fasta_resample.read()
                if event_resample is None:
                    break

                # get the input values from the user
                no_seqs = values_resample['seqs']
                perc = values_resample['seqs']
                iterations = values_resample['iter']

######################################################################################################################
# RUN THE RESAMPLING ANALYSIS
######################################################################################################################

                if event_resample == 'Run':

                    if values_resample['set_seed']:
                        random.seed(1234)

                    if values_resample['inc_outgroups'] and len(outgroup_display) == 0:
                        sg.popup_ok('Outgroups?', 'You have indicated that outgroups should be included, '
                                                  'but you have not uploaded a Fasta file for them.')

                    elif len(no_seqs) == 0 or len(iterations) == 0:
                        sg.popup_ok('Iterations?', 'You have not specified values for the number of '
                                                   'sequences and/or iterations.')

                    elif int(no_seqs) <= 0 or int(iterations) <= 0:
                        sg.popup_ok('Negative or Zero', 'You have inserted negative numbers and/or zeroes.')

                    else:
                        save_dir = values_resample['Browse']
                        os.chdir(save_dir)
                        # This if statement checks if the Iterations folder already exists, and then overwrites it
                        # if it does
                        if os.path.exists(save_dir + '\\' + 'Iterations' + '_' + perc + '/'):
                            shutil.rmtree(save_dir + '\\' + 'Iterations' + '_' + perc + '/')

                        os.mkdir('Iterations' + '_' + perc + '/')

                        f1 = ()
                        counter = 0

                        no_seqs = (int(values_resample['seqs']) / 100) * n
                        no_seqs = round(no_seqs, 0)

                        while counter < int(iterations):

                            r_shuffle = random.sample(content2, k=int(no_seqs))

                            for line in r_shuffle:
                                # splitting the line here takes the fasta file back to its original form,
                                # with the name on one line, and the sequence below it
                                y = line.split()
                                f1 = open('Iterations' + '_' + perc + '/iteration' + str(counter + 1) + ".fas", "a")
                                for element in y:
                                    f1.write(element + "\n")
                            # if outgroups is selected, append them to the end of each file
                            if values_resample['inc_outgroups']:
                                outgroups = outgroup_display.splitlines()
                                for lines in outgroups:
                                    f1.write(lines + "\n")
                                # noinspection PyUnboundLocalVariable
                                outgroup_fasta_file.close()
                            counter += 1
                            f1.close()

                            window_fasta_resample['progbar_fasta'].update_bar(counter, max=int(iterations))

                        sg.popup_ok(str(iterations) + ' iteration files successfully written to ' + save_dir +
                                    '/Iterations_' + perc + '. Each iteration file contains ' + str(no_seqs)
                                    + ' sequences.')

######################################################################################################################
# CONVERT FASTA TO PHYLIP FOR RAXML
######################################################################################################################

    if event == 'Convert Fasta to Phylip':
        convert_folder = sg.popup_get_folder('fasta files to convert', no_window=True)

        if len(convert_folder) > 0:
            full_filenames_raxml = []
            short_filename_raxml = []

            for file in os.listdir(convert_folder):
                name = os.fsdecode(file)
                if name.endswith((".fas", ".fasta", ".FASTA")):
                    full_filenames_raxml.append(convert_folder + "\\" + name)
                    short_filename_raxml.append(name)

            if len(full_filenames_raxml) == 0:
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

                for f in full_filenames_raxml:
                    sg.one_line_progress_meter('Progress', current_value=i, max_value=len(full_filenames_raxml),
                                               grab_anywhere=True, key='conversion_progress')
                    if f.endswith(".fas"):
                        AlignIO.convert(f, "fasta", "Phylip_conversions/" + short_filename_raxml[i].replace(".fas", "")
                                        + ".phy", "phylip", "DNA")
                    elif f.endswith(".fasta"):
                        AlignIO.convert(f, "fasta",
                                        "Phylip_conversions/" + short_filename_raxml[i].replace(".fasta", "") + ".phy",
                                        "phylip", "DNA")
                    elif f.endswith(".FASTA"):
                        AlignIO.convert(f, "fasta",
                                        "Phylip_conversions/" + short_filename_raxml[i].replace(".FASTA", "") + ".phy",
                                        "phylip", "DNA")
                    i += 1

                sg.one_line_progress_meter_cancel(key='conversion_progress')
                sg.popup_ok(str(len(full_filenames_raxml)) + " files successfully converted from .fas to .phy. See the "
                                                       "\'Phylip_conversions\' folder in the " + convert_folder
                            + ' directory.')

            except ValueError:
                # sequence names must not be more than 10 chars, otherwise the rest is cut off, and then there may be
                # repeated names
                sg.popup_ok('File Issue', 'Unequal sequence lengths or repeated names. '
                                          'Please check your file and re-upload.')

######################################################################################################################
# end of try-catch block
######################################################################################################################

######################################################################################################################
# RUN RAxML for a single input fasta file
######################################################################################################################

    if event == 'Run RAxML (single Phylip file)':
        layout_raxml_single = [
            [sg.Text('Select your Fasta file:')],
            [sg.Input(), sg.FileBrowse(key='raxml_infile')],
            [sg.Text('Model:')],
            [sg.InputCombo(('GTRCAT', 'GTRGAMMA', 'GTRMIX', 'GTRCAT_GAMMA', 'GTRGAMMAI'),
                           size=(20, 1), key='raxml_model_single')],
            [sg.Text('Number of bootstrap iterations (please insert a value greater than zero):')],
            [sg.InputText(key='bootstraps_raxml_single', size=(10, 10))],
            [sg.Checkbox('Set seed?', key='raxml_seed_single', default=False, tooltip='Select this if you want '
                            'the results to be reproducible if you run this again.')],
            [sg.Submit('Run', key='run_raxml_single')]
            ]

        # open a new window with settings
        window_raxml_single = sg.Window('RAxML', layout_raxml_single, size=(500, 300))

        while True:
            event, values = window_raxml_single.read()
            if event is None:
                break

            raxml_input_file = values['raxml_infile']
            boots_raxml = values['bootstraps_raxml_single']
            model = values['raxml_model_single']

            if event == 'run_raxml_single':

                try:
                    boots_raxml_int = int(boots_raxml)
                except ValueError:
                    sg.popup_ok('Error', 'Please input an integer for bootstrap number.')
                    continue

                if len(raxml_input_file) == 0:
                    sg.popup_ok('Select File', 'Please select a Phylip file.')
                elif len(model) == 0:
                    sg.popup_ok('Model Selection', 'Please select an evolutionary model from the dropdown menu.')
                elif len(boots_raxml) == 0:
                    sg.popup_ok('Bootstrap Value', 'Number of bootstraps missing. Please check.')
                elif int(boots_raxml) < 0:
                    sg.popup_ok('Negative Value', 'Please insert a positive bootstrap value.')

                else:

                    if os.path.exists(os.path.dirname(raxml_input_file) + '\\' + 'RAxML_output_single/'):
                        shutil.rmtree(os.path.dirname(raxml_input_file) + '\\' + 'RAxML_output_single/')

                    os.mkdir(os.path.dirname(raxml_input_file) + '\\' + 'RAxML_output_single/')

                    try:
                        # specify the name of the output tree file:
                        raxml_outfile = os.path.basename(raxml_input_file).replace(".phy", "")

                        os.chdir(raxmlHPC_path)

                        # if the user sets a seed:
                        if values['raxml_seed_single']:
                            cmd_raxml = RaxmlCommandline(sequences=raxml_input_file, name=raxml_outfile,
                                                         model=model,
                                                         num_replicates=int(boots_raxml),
                                                         arsimony_seed=int(values['raxml_seed_single']),
                                                         working_dir=os.path.dirname(raxml_input_file) +
                                                                     '\\' + 'RAxML_output_single/')
                        else:
                            cmd_raxml = RaxmlCommandline(sequences=raxml_input_file, name=raxml_outfile,
                                                         model=model,
                                                         num_replicates=int(boots_raxml),
                                                         working_dir=os.path.dirname(raxml_input_file) + '\\' + 'RAxML_output_single/')

                        cmd_raxml()

                        #window_raxml_single['progbar_raxml'].update_bar(j, max=len(full_filenames_raxml))

                        sg.popup_ok('Success! ML tree successfully created using ' + str(boots_raxml) +
                                    ' bootstrap repeats. These have been written to the ' + raxml_input_file +
                                    '\\' + 'RAxML_output folder.')

                    except:
                        sg.popup_ok('RAxML encountered an error. This could be due to (1) your input '
                                    'files having too few sequences, or very short sequence lengths, '
                                    'or (2) an issue with the model you have chosen.')

######################################################################################################################
# RUN RAxML for multiple files in a folder
######################################################################################################################

    if event == 'Run RAxML (multiple Phylip files)':
        layout_raxml = [
            [sg.Text('Select the folder containing your resampled Fasta files:')],
            [sg.Input(), sg.FolderBrowse(key='raxml_infolder')],
            [sg.Text('Model:')],
            [sg.InputCombo(('GTRCAT', 'GTRGAMMA', 'GTRMIX', 'GTRCAT_GAMMA', 'GTRGAMMAI'),
                           size=(20, 1), key='raxml_model')],
            [sg.Text('Number of bootstrap iterations (please insert a value greater than zero):')],
            [sg.InputText(key='bootstraps_raxml', size=(10, 10))],
            [sg.Checkbox('Set seed?', key='raxml_seed', default=False, tooltip='Select this if you want '
                            'the results to be reproducible if you run this again.')],
            [sg.Submit('Run', key='run_raxml')],
            [sg.ProgressBar(max_value=len(full_filenames_raxml), orientation='h', size=(20, 20), key='progbar_raxml')]
        ]
        # open a new window with settings
        window_raxml = sg.Window('RAxML', layout_raxml, size=(500, 300))

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
                    full_filenames_raxml = []
                    short_filenames_raxml = []

                    if os.path.exists(raxml_input_folder + '\\' + 'RAxML_output/'):
                        shutil.rmtree(raxml_input_folder + '\\' + 'RAxML_output/')

                    os.mkdir(raxml_input_folder + '\\' + 'RAxML_output/')

                    for file in os.listdir(raxml_input_folder):
                        name = os.fsdecode(file)
                        if name.endswith(".phy"):
                            full_filenames_raxml.append(raxml_input_folder + "\\" + name)
                            short_filenames_raxml.append(name)

                    if len(full_filenames_raxml) == 0:
                        sg.popup_ok("There are no .phy files in the selected folder.")

                    else:
                        try:
                            j = 0
                            # Loop through each file in the desired folder
                            for file in full_filenames_raxml:
                                raxml_input = file
                                # specify the name of the output tree file:
                                raxml_outfile = short_filenames_raxml[j].replace(".phy", "")

                                os.chdir(raxmlHPC_path)

                                # if the user sets a seed:
                                if values['raxml_seed']:
                                    cmd_raxml = RaxmlCommandline(sequences=raxml_input, name=raxml_outfile, model=model,
                                                                 num_replicates=int(boots_raxml),
                                                                 parsimony_seed=int(values['raxml_seed']),
                                                                 working_dir=raxml_input_folder +
                                                                             '\\' + 'RAxML_output/')
                                else:
                                    cmd_raxml = RaxmlCommandline(sequences=raxml_input, name=raxml_outfile, model=model,
                                                                 num_replicates=int(boots_raxml),
                                                                 working_dir=raxml_input_folder + '\\' + 'RAxML_output/')

                                cmd_raxml()
                                j += 1
                                window_raxml['progbar_raxml'].update_bar(j, max=len(full_filenames_raxml))

                            sg.popup_ok('Success!', str(len(full_filenames_raxml)) +
                                        ' ML trees successfully created using ' + str(boots_raxml) +
                                        ' bootstrap repeats. This has been written to the ' + raxml_input_folder +
                                        '\\' + 'RAxML_output folder.')

                        except:
                            sg.popup_ok('RAxML encountered an error. This could be due to (1) your input '
                                        'files having too few sequences, or very short sequence lengths, '
                                        'or (2) an issue with the model you have chosen.')

######################################################################################################################
# FASTTREE for multiple fasta files in a folder
######################################################################################################################

    if event == 'Run FastTree (multiple Fasta files)':
        layout_fasttree = [
            [sg.Text('Select the folder containing your resampled Fasta files:')],
            [sg.Input(), sg.FolderBrowse(key='fasttree_infolder')],
            [sg.Text('Number of bootstrap iterations:')],
            [sg.InputText(key='bootstraps_fasttree', size=(10, 10))],
            [sg.Checkbox('GTR', key='gtr', tooltip='Apply the generalized time-reversible model \ninstead of the '
                                                   'default Jukes-Cantor model \n(nucleotide data only)'),
             sg.Checkbox('Gamma', key='gamma', tooltip='Apply the discrete gamma model'),
             sg.Checkbox('Nucleotide', default=True, key='nt', tooltip='Deselect for protein data'),
             sg.Checkbox('Seed?', default=False, key='fasttree_seed', tooltip='Select this if you want '
                            'the resampling results to be reproducible if you run this again.')],
            [sg.Submit('Run', key='run_fasttree')],
            [sg.ProgressBar(max_value=len(full_filenames_fasttree), orientation='h', size=(20, 20),
                            key='progbar_fasttree')]
        ]
        # open a new window with settings
        window_fasttree = sg.Window('FastTree', layout_fasttree, size=(500, 250))

        while True:
            event_fasttree, values_fasttree = window_fasttree.read()
            if event_fasttree is None:
                break

            if event_fasttree == 'run_fasttree':

                try:
                    boots_fasttree = values_fasttree['bootstraps_fasttree']
                    boots_fasttree_int = int(boots_fasttree)
                except ValueError:
                    sg.popup_ok('Error', 'Bootstrap value must be an integer.')
                    continue

                fasttree_input_folder = values_fasttree['fasttree_infolder']

                if len(fasttree_input_folder) == 0:
                    sg.popup_ok('Select Folder', 'Please select a folder containing your Fasta files.')
                elif len(boots_fasttree) == 0:
                    sg.popup_ok('Bootstrap Value', 'Please input the number of bootstraps you wish to run.')
                elif int(boots_fasttree) < 0:
                    sg.popup_ok('Negative Value', 'Please insert a positive bootstrap value.')
                else:
                    full_filenames_fasttree = []
                    short_filenames_fasttree = []

                    for file in os.listdir(fasttree_input_folder):
                        name = os.fsdecode(file)
                        if name.endswith(".fas") or name.endswith(".fasta"):
                            full_filenames_fasttree.append(fasttree_input_folder + "\\" + name)
                            short_filenames_fasttree.append(name)

                    if len(full_filenames_fasttree) == 0:
                        sg.popup_ok("There are no Fasta files in the selected folder.")

                    else:
                        # create an output folder to save trees to
                        if os.path.exists(fasttree_input_folder + '\\' + 'FastTree_output/'):
                            shutil.rmtree(fasttree_input_folder + '\\' + 'FastTree_output/')

                        os.mkdir(fasttree_input_folder + '\\' + 'FastTree_output/')

                    try:
                        q = 0
                        # Loop through each file in the desired folder
                        for file in full_filenames_fasttree:
                            fasttree_fastafile = file
                            # specify the location of the output tree file:
                            if file.endswith('.fas'):
                                fasttree_outfile = fasttree_input_folder + '\\' + 'FastTree_output/' + \
                                                    short_filenames_fasttree[q].replace('.fas', '.tre')
                            elif file.endswith('.FASTA'):
                                fasttree_outfile = fasttree_input_folder + '\\' + 'FastTree_output/' + \
                                                    short_filenames_fasttree[q].replace('.FASTA', '.tre')
                            elif file.endswith('.fasta'):
                                fasttree_outfile = fasttree_input_folder + '\\' + 'FastTree_output/' + \
                                                    short_filenames_fasttree[q].replace('.fasta', '.tre')
                            # Run FastTree:
                            # if the user wants a seed:
                            if values_fasttree['fasttree_seed']:
                                fasttree_seed = 1234
                                cmd_fasttree = FastTreeCommandline(fasttree_exe, nt=values_fasttree['nt'], gtr=values_fasttree['gtr'],
                                                                   gamma=values_fasttree['gamma'], input=fasttree_fastafile,
                                                                   out=fasttree_outfile, boot=int(boots_fasttree),
                                                                   seed=int(values_fasttree['fasttree_seed']))
                            # if no seed required:
                            else:
                                cmd_fasttree = FastTreeCommandline(fasttree_exe, nt=values_fasttree['nt'], gtr=values_fasttree['gtr'],
                                                                   gamma=values_fasttree['gamma'], input=fasttree_fastafile,
                                                                   out=fasttree_outfile, boot=int(boots_fasttree))
                            cmd_fasttree()
                            q += 1
                            window_fasttree['progbar_fasttree'].update_bar(q, max=len(full_filenames_fasttree))

                        sg.popup_ok('Success!', 'ML trees successfully created with ' + str(boots_fasttree) +
                                    ' bootstrap repeats. These have been written to the ' + fasttree_input_folder +
                                    '\\' + 'FastTree_output folder.')

                    except:
                        sg.popup_ok('FastTree encountered an error.')

######################################################################################################################
# FASTTREE for a single fasta file
######################################################################################################################

    if event == 'Run FastTree (single Fasta file)':
        layout_fasttree_single = [
            [sg.Text('Select the file containing your Fasta file:')],
            [sg.Input(), sg.FileBrowse(key='fasttree_infile')],
            [sg.Text('Number of bootstrap iterations:')],
            [sg.InputText(key='bootstraps_fasttree_single', size=(10, 10))],
            [sg.Checkbox('GTR', key='gtr', tooltip='Apply the generalized time-reversible model \ninstead of the '
                                                   'default Jukes-Cantor model \n(nucleotide data only)'),
            sg.Checkbox('Gamma', key='gamma', tooltip='Apply the discrete gamma model'),
            sg.Checkbox('Nucleotide', default=True, key='nt', tooltip='Deselect for protein data'),
            sg.Checkbox('Seed?', default=False, key='fasttree_seed_single')],
            [sg.Submit('Run', key='run_fasttree_single')],
            ]
        # open a new window with settings
        window_fasttree_single = sg.Window('FastTree', layout_fasttree_single, size=(500, 250))

        while True:
            event, values = window_fasttree_single.read()
            if event is None:
                break

            fasttree_input_file = values['fasttree_infile']
            boots_fasttree = values['bootstraps_fasttree_single']

            if event == 'run_fasttree_single':

                try:
                    boots_fasttree_int = int(boots_fasttree)
                except ValueError:
                    sg.popup_ok('Error', 'Bootstrap value must be an integer.')
                    continue

                if len(fasttree_input_file) == 0:
                    sg.popup_ok('Select File', 'Please select a Fasta file.')
                elif len(boots_fasttree) == 0:
                    sg.popup_ok('Bootstrap Value', 'Please input the number of bootstraps you wish to run.')
                elif int(boots_fasttree) < 0:
                    sg.popup_ok('Negative Value', 'Please insert a positive bootstrap value.')
                else:
                    # create an output folder to save trees to
                    if os.path.exists(os.path.dirname(fasttree_input_file) + '\\' + 'FastTree_single/'):
                        shutil.rmtree(os.path.dirname(fasttree_input_file) + '\\' + 'FastTree_single/')

                    os.mkdir(os.path.dirname(fasttree_input_file) + '\\' + 'FastTree_single/')

                try:
                    # specify the location of the output tree file:
                    if fasttree_input_file.endswith('.fas'):
                        fasttree_outfile = os.path.dirname(fasttree_input_file) + '\\' + 'FastTree_single/' + \
                                                       os.path.basename(fasttree_input_file).replace('.fas', '.tre')
                    elif fasttree_input_file.endswith('.FASTA'):
                        fasttree_outfile = os.path.dirname(fasttree_input_file) + '\\' + 'FastTree_single/' + \
                                                       os.path.basename(fasttree_input_file).replace('.FASTA', '.tre')
                    elif fasttree_input_file.endswith('.fasta'):
                        fasttree_outfile = os.path.dirname(fasttree_input_file) + '\\' + 'FastTree_single/' + \
                                                       os.path.basename(fasttree_input_file).replace('.fasta', '.tre')
                    # Run FastTree:
                    # if the user wants a seed:
                    if values['fasttree_seed_single']:
                        fasttree_seed = 1234
                        cmd_fasttree = FastTreeCommandline(fasttree_exe, nt=values_fasttree['nt'],
                                                           gtr=values['gtr'],
                                                           gamma=values['gamma'],
                                                           input=fasttree_input_file,
                                                           out=fasttree_outfile, boot=int(boots_fasttree),
                                                           seed=int(values['fasttree_seed_single']))
                    # if no seed required:
                    else:
                        cmd_fasttree = FastTreeCommandline(fasttree_exe, nt=values['nt'],
                                                           gtr=values['gtr'],
                                                           gamma=values['gamma'],
                                                           input=fasttree_input_file,
                                                           out=fasttree_outfile, boot=int(boots_fasttree))

                    cmd_fasttree()

                    sg.popup_ok('Success! ML tree successfully created with ' + str(boots_fasttree) +
                                ' bootstrap repeats. This has been written to the ' + os.path.basename(fasttree_input_file) +
                                '\\' + 'FastTree_single folder.')

                except:
                    sg.popup_ok('FastTree encountered an error.')

######################################################################################################################
# CONTACT DETAILS WINDOW
######################################################################################################################

    if event == 'Contact':
        sg.popup_scrolled('Report bugs or make suggestions',  'Clarke van Steenderen \nemail: vsteenderen@gmail.com or '
                                                        'g14v1511@campus.ru.ac.za')

######################################################################################################################
# CITATION DETAILS WINDOW
######################################################################################################################

    if event == 'Citation':
        sg.popup_scrolled('Please cite this program as follows: \nvan Steenderen, CJM. SPEDE-SAMPLER: '
                            'assess sampling effects on species delimitation, version 1.1, 2020.')

window_main.close()
