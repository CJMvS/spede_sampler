import PySimpleGUI as sg
import os
import random
import shutil
from Bio import AlignIO
from Bio.Phylo.Applications import RaxmlCommandline

sg.theme('DarkBlack 1')

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
global boots
boots = 0

######################################################################################################################
# MENU OPTIONS
######################################################################################################################

menu_def = [
           ['Resample Fasta files', ['Open Ingroup Fasta File', 'Print Ingroup Fasta File to Screen',
            'Open Outgroup Fasta File', 'Print Outgroup Fasta File to Screen', 'Nucleotide Composition',
                                     'Run Resampling', 'Fasta Resampling Help', 'Exit']],
           ['RAxML', ['Convert Fasta to Phylip', 'Run RAxML', 'RAxML Help']],
           ['FastTree', ['Select Folder with Resampled Fasta Files', 'Run FastTree', 'FastTree Help']],
           ['GMYC', ['Open Folder with ML phylogenies', 'Run GMYC', 'Help']],
           ['About', ['Citation', 'Contact']], ]

######################################################################################################################
# MAIN WINDOW LAYOUT
######################################################################################################################

layout_main = [
    [sg.Menu(menu_def, tearoff=False)],
    [sg.Image('spede_sampler.png')]
]

######################################################################################################################
# CREATE MAIN WINDOW
######################################################################################################################

window_main = sg.Window('SPEDE-RESAMPLER', layout_main, size=(560, 350), icon='bug.ico')

while True:
    event, values = window_main.read()

    if event is None:
        break

######################################################################################################################
# HELP WITH FASTA FILE INPUT FOR RESAMPLING
######################################################################################################################

    if event == 'Fasta Resampling Help':
        help_fasta_resampling = open('help_fasta_resampling.txt', 'r')
        help_fasta_resampling_content = help_fasta_resampling.read()
        help_fasta_resampling.close()
        sg.popup_scrolled(help_fasta_resampling_content, title='Help', size=(100, 25))

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
            sg.popup_notify('Please upload an ingroup Fasta file to display',
                            location=(600, 400), display_duration_in_ms=4000,
                            fade_in_duration=10, alpha=1)
        else:
            # noinspection PyUnboundLocalVariable
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
            sg.popup_notify('Please upload an outgroup Fasta file to display',
                            location=(600, 400), display_duration_in_ms=4000,
                            fade_in_duration=10, alpha=1)
        else:
            # noinspection PyUnboundLocalVariable
            sg.popup_scrolled(outgroup_display, title=outgroup_fasta_filename, size=(100, 100))

######################################################################################################################
# NUCLEOTIDE COMPOSITION
######################################################################################################################

    if event == 'Nucleotide Composition':

        if len(ingroup_display) == 0:
            sg.popup_notify('You have not uploaded an ingroup Fasta file', location=(600, 400),
                            display_duration_in_ms=5000, fade_in_duration=10, alpha=1)
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
                                line[i] == "B" or line[i] == "D" or line[i] == "H" or line[i] == "V":
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
            sg.popup_notify('You have not uploaded an ingroup Fasta file', location=(600, 400),
                            display_duration_in_ms=5000, fade_in_duration=10, alpha=1)

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
                        sg.popup_notify('You have indicated that outgroups should be included, but you have not '
                                        'uploaded a Fasta file for them.',
                                        location=(550, 400), display_duration_in_ms=5000,
                                        fade_in_duration=10, alpha=1)

                    elif len(no_seqs) == 0 or len(iterations) == 0:
                        sg.popup_notify('You have not specified values for the number of sequences and iterations',
                                        location=(550, 400), display_duration_in_ms=5000,
                                        fade_in_duration=10, alpha=1)

                    elif int(no_seqs) <= 0 or int(iterations) <= 0:
                        sg.popup_notify('You have inserted negative numbers', location=(650, 400),
                                        display_duration_in_ms=5000, fade_in_duration=10, alpha=1)

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
                sg.popup_notify("There are no FASTA files in the selected folder.", location=(650, 400),
                                display_duration_in_ms=6000, fade_in_duration=10, alpha=1)

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
                sg.popup_notify('Unequal sequence lengths or repeated names. Please check your file and re-upload.',
                                location=(650, 400), display_duration_in_ms=6000, fade_in_duration=10, alpha=1)

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
            [sg.Input(), sg.FolderBrowse('Input Folder')],
            [sg.Text('Model:')],
            [sg.InputCombo(('GTRCAT', 'GTRCAT_GAMMAI', 'GTRCAT_GAMMA', 'GTRGAMMAI', 'GTRGAMMA', 'GTRMIX', 'GTRMIXI'),
                           size=(20, 1), key='raxml_model')],
            [sg.Text('Number of bootstrap iterations:')],
            [sg.InputText(key='bootstraps', size=(10, 10))],
            [sg.Submit('Run Analysis')],
            [sg.ProgressBar(10, orientation='h', size=(20, 20), key='progbar_raxml')]
        ]
        # open a new window with settings
        window_raxml = sg.Window('RAxML', layout_raxml, size=(500, 250))

        while True:
            event, values = window_raxml.read()
            if event is None:
                break

            raxml_input_folder = values['Input Folder']
            boots = values['bootstraps']
            model = values['raxml_model']

            if event == 'Run Analysis':

                if len(raxml_input_folder) == 0:
                    sg.popup_notify('Please select a folder containing your Phylip files.',
                                    location=(650, 400), display_duration_in_ms=6000, fade_in_duration=10, alpha=1)
                elif len(boots) == 0:
                    sg.popup_notify('Please input the number of bootstraps you wish to run.',
                                    location=(650, 400), display_duration_in_ms=6000, fade_in_duration=10, alpha=1)
                elif len(model) == 0:
                    sg.popup_notify('Please select an evolutionary model from the dropdown menu.',
                                    location=(650, 400), display_duration_in_ms=6000, fade_in_duration=10, alpha=1)

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
                        sg.popup_notify("There are no .phy files in the selected folder.")
                    else:
                        j = 0
                        # Loop through each file in the desired folder
                        for file in full_filenames:
                            raxml_input = file
                            # specify the name of the output tree file:
                            raxml_outfile = short_filenames[j].replace(".phy", "")
                            cmd_raxml = RaxmlCommandline(sequences=raxml_input, name=raxml_outfile, model=model,
                                                         num_replicates=boots,
                                                         working_dir=raxml_input_folder + '\\' + 'RAxML_output/')
                            cmd_raxml()
                            j += 1
                            window_raxml['progbar_raxml'].update_bar(j + 1)

                        sg.popup_notify(str(boots) + ' ML trees successfully created. These have been written to the '
                                        + raxml_input_folder + '\\' + 'RAxML_output folder.', location=(650, 400),
                                        display_duration_in_ms=6000, fade_in_duration=10, alpha=1)

######################################################################################################################
# EXIT THE PROGRAM
######################################################################################################################

    if event == 'Exit':
        window_main.close()

window_main.close()
