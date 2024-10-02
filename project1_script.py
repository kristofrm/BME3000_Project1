# -*- coding: utf-8 -*-

"""
BME3000 Project 1
Visualizing ECG Data of Normal and Arrhythmic Events
Kristof Rohaly-Medved & Sabrina Goslin
Due 10/03/24
"""

#Importing necessary packages
import numpy as np
from matplotlib import pyplot as plt
import project1_module as p1m

#Naming the file
filename = 'ecg_e0103_half1.npz'

#Call the function
ecg_voltage, fs, label_samples, label_symbols, subject_id, electrode, units = p1m.load_data(filename)
dt = 1/fs #seconds per sample, 1 sample = 1/250th of a second (250 samples = 1 second)
 
#Create a time array to hold ecg data
ecg_time = np.arange(0,3600,dt)
samples_per_second = 250

#%% Part 2: Plotting the raw signal
 
#Creating a figure
plt.figure(1, clear=True)

#Calling the function and displaying the plot
p1m.plot_raw_data(ecg_voltage, ecg_time, 'V', f'ECG Voltage vs. Time\nSubject ID: {subject_id} for electrode {electrode}')

#%% Part 3: Plotting events

#Call the function, set x limits to display two normal events and one arrhythmic event
p1m.plot_events(label_samples, label_symbols, ecg_time, ecg_voltage)
plt.xlim(1890,1900)

#%% Part 4: Extracting Trials

#Use np.unique() to get unique label symbols
unique_labels = np.unique(label_symbols)

#Split labels into N and V using indexes of unique_labels array
event_samples_n = label_samples[label_symbols == unique_labels[0]]
event_samples_v = label_samples[label_symbols == unique_labels[1]]
 
#Create an array filled with zeros as long as the amount of normal events
n_start_indices = np.zeros(len(event_samples_n), dtype=int)
 
#Loop through each value to get starting indices for each event (event_index - 250/2)
for i, event_index in enumerate(event_samples_n): #Got idea to use enumerate function from ChatGPT, modified on our own.

    n_start_indices[i] = max(event_index - samples_per_second/2, 0) # Got idea to use max function from ChatGPT, implemented indepenedently.
    
#Create an array filled with zeros as long as the amount of arrhythmic events
v_start_indices = np.zeros(len(event_samples_v), dtype=int)
 
#Loop through each value to get starting indices for each event (event_index - 250/2)
for i, event_index in enumerate(event_samples_v):
    v_start_indices[i] = max(event_index - samples_per_second/2, 0)
    
#Call function for n and v event
n_trials = p1m.extract_trials(ecg_voltage, n_start_indices, samples_per_second)
v_trials = p1m.extract_trials(ecg_voltage, v_start_indices, samples_per_second)

#Using an if statement, verify they are the expected size and are full of values. Print a message if they are
#Check for both event_samples_n and event_samples_v in a single if statement
if len(event_samples_n) == len(n_trials[:,0]) and len(event_samples_v) == len(v_trials[:,0]):
    n_valid = event_samples_n.size > 0 and np.all(np.isfinite(event_samples_n))
    v_valid = event_samples_v.size > 0 and np.all(np.isfinite(event_samples_v))
    if n_valid and v_valid:
        print(f'There are {len(event_samples_n)} N events and {len(n_trials[:, 0])} rows in n_trials.')
        print(f'There are {len(event_samples_v)} V events and {len(v_trials[:, 0])} rows in v_trials.')
        print('The trial arrays match up correctly.')
    else:
        if not n_valid:
            print('event_samples_n is empty or contains invalid values.')
        if not v_valid:
            print('event_samples_v is empty or contains invalid values.')
else:
    print('The lengths of event_samples_n and event_samples_v do not match the respective trial arrays.')
print()

#Specific time of event is arbitrary. Used 0-1 second interval to display both events together.
trial_times_plot = np.arange(-0.5,0.5,dt)

#Creating a second plot of data from trials N and V
plt.figure(2, clear=True)

#Plot the first normal and abnormal events
plt.plot(trial_times_plot, n_trials[1, :], label='Normal Event')
plt.plot(trial_times_plot, v_trials[1, :], label='Abnormal Event')  

#Set labels and title
plt.xlabel('Event Time (seconds)')
plt.ylabel(f'Event Voltage ({units})')
plt.title('Event Time vs. Voltage for Normal and Abnormal Events')
plt.grid(True)
plt.legend()

#%% Part 5: Plot Trial Means

#Set duration of each 'trial'
trial_duration_seconds = 1 # each event goes from -0.5s before to 0.5s after event

symbols, trial_time, mean_trial_signal = p1m.plot_mean_and_std_trials(ecg_voltage, label_samples, label_symbols, trial_duration_seconds, fs, units, f'Mean and STD of ECG Voltage vs. Time\n Comparing Normal and Abnormal Cardiac Events\nSubject ID: {subject_id} for electrode {electrode}')

#%% Part 6: Save arrays and plots
 
out_filename = f'ecg_means_{subject_id}.npz'
    
# Call function to save variables to npz file
p1m.save_means(symbols, trial_time, mean_trial_signal, out_filename)

# Load in previously saved npz file and extract variables
infile = np.load(out_filename)
symbols_loaded = infile['symbols']
trial_time_loaded = infile['trial_time']
mean_trial_signal_loaded = infile['mean_trial_signal']

# Compare loaded arrays to existing arrays
if np.array_equal(symbols_loaded, symbols):
    print('The loaded symbols array matches the existing symbols array.')
else:
    print('The loaded symbols array does not match the exsting symbols array.')
if np.array_equal(trial_time_loaded, trial_time):
    print('The loaded trial time array matches the existing trial time array.')
else:
    print('The loaded trial time array does not matche the existing trial time array.')
if np.array_equal(mean_trial_signal_loaded, mean_trial_signal):
    print('The loaded mean trial signal array matches the existing mean trial signal array.')
else:
    print('The loaded mean trial signal array does not match the existing mean trial signal array.')
    
# Saving the figures
plt.figure(1)
plt.savefig(f'ecg_voltage_over_time_events_labeled_subject_{subject_id}.png')
plt.figure(2)
plt.savefig(f'ecg_voltage_over_time_for_single_n&v_event_subject_{subject_id}.png')
plt.figure(3)
plt.savefig(f'ecg_voltage_over_time_mean_and_std_for_n&v_events_subject_{subject_id}.png')