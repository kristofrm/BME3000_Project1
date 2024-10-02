#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 12:08:07 2024

@author: kristofrm
"""

#Importing necessary packages
import numpy as np
from matplotlib import pyplot as plt

#%% Part 1: Load data
    
#Creating a function to load the data and extract necessary attributes
def load_data(input_file):
    """
    Parameters
    ----------
    input_file : file
        A file containing multiple types of information.
 
    Returns
    -------
    ecg_voltage : array of floats
        Array of actual ECG values at each sample.
    fs : int
        The frequency at which samples are taken in a second (250).
    label_samples : array of ints
        Indexes at which events occur.
    label_symbols : array of strs
        An N or V, representing normal or arrhythmetic events in the whole data set.
    subject_id : str
        Gives the ID of the person being monitored.
    electrode : str
        Gives the type of electrode used.
    units : str
        Provides units of volts, given as mV.
 
    """
    # Load data file and print files
    data = np.load(input_file)
    print(data.files)
    print()
    
    # Extract relevant data files
    ecg_voltage = data['ecg_voltage']
    fs = data['fs']
    label_samples = data['label_samples']
    label_symbols = data['label_symbols']
    subject_id = data['subject_id']
    electrode = data['electrode']
    units = data['units']
    
    return ecg_voltage, fs, label_samples, label_symbols, subject_id, electrode, units

#%% Part 2: Plotting the raw signal

#Writing a function to plot the raw data
def plot_raw_data(signal_voltage, signal_time, units='V', title=''):
    """
    Parameters
    ----------
    signal_voltage : array of floats
        Array of size ecg_voltage, gives voltage values at given samples.
    signal_time : array of floats
        Array of size ecg_voltage, puts time in seconds on the x-axis.
    units : str, optional
        Gives the units of ECG values. The default is 'V'.
    title : str, optional
        Provides a title for the plot. The default is ''.
 
    Returns
    -------
    None.
 
    """
    plt.plot(signal_time, signal_voltage, label=f'Signal ({units})', color='b')
    
    #Annotating the plot
    plt.xlabel('Time (seconds)')
    plt.ylabel(f'Voltage ({units})')
    plt.title(title)
    plt.grid(True)
    plt.legend()
    
    
#%% Part 3: Plotting events
    
#Define a function to plot the events
def plot_events(label_samples, label_symbols, signal_time, signal_voltage):
    """
    Parameters
    ----------
    label_samples : array of ints
        1D array of 3661 indexes at which events occur.
    label_symbols : array of strings
        An array size of label_samples filled with N or V, representing normal or arrhythmetic events in the whole data set.
    signal_time : float
        Time in seconds.
    signal_voltage : array of floats
        Array size of ecg_data, gives voltage values at given samples.
 
    Returns
    -------
    None.
 
    """
    
    #Use np.unique() to get unique label symbols
    unique_labels = np.unique(label_symbols)
    
    #Use a for loop to split labels in 'N' and 'V'
    for label in unique_labels:
        event_samples_label = label_samples[label_symbols == label]
        signal_time_label = signal_time[event_samples_label]
        signal_voltage_label = signal_voltage[event_samples_label]
        
        #Plot each event as a dot
        plt.scatter(signal_time_label, signal_voltage_label, label=f'Event {label}', zorder=2)
        
    #Give the plot a legend in the lower right corner
    plt.legend(loc='lower right')
    
#%% Part 4: Extracting Trials
    
#Define a function to extract trials
#Create a function to extract trials
#note: trial_sample_count = # samples each trial should have (# columns)
def extract_trials(signal_voltage, trial_start_samples, trial_sample_count):
    """
    Parameters
    ----------
    signal_voltage : array
        Array of size ecg_data that gives signal voltages at given samples.
    trial_start_samples : TYPE
        DESCRIPTION.
    trial_sample_count : int = 250
        Number of samples per 1 second interval (equal to sampling frequency). Specifies number of columns for trials array
    Returns
    -------
    trials : array of floats
        2D array of ecg voltage values. Rows = number of events, columns = ecg voltage sample
 
    """
    
    trials = np.zeros((len(trial_start_samples), trial_sample_count), dtype=float)
    # Loop through each event starting index in trial_start_samples and add 
    # ecg starting, event, and end data to trials array based off trial_start_samples index
    for i, event_start_index in enumerate(trial_start_samples):
        if not event_start_index + trial_sample_count > len(signal_voltage):
            trials[i, :] = signal_voltage[event_start_index:event_start_index + trial_sample_count]
            
    return trials

#%% Part 5: Plot Trial Means

#Create a function that plots means and standard deviations of ab/normal events
def plot_mean_and_std_trials(signal_voltage, label_samples, label_symbols, trial_duration_seconds, fs, units, title):
    """
    Parameters
    ----------
    signal_voltage : array
        An aarray of size ecg_data filled with the value of signal voltage at given sample.
    label_samples : array
        An array of the size of ecg_data filled with indices at which events occur.
    label_symbols : array
        An array the size of ecg_data filled with N or V representing normal or arrhythmic events.
    trial_duration_seconds : int
        Value of how long the trial lasts (1 second).
    fs : int
        The frequency at which samples are taken per second (250).
    units : str
        Units of ECG data (mV).
    title : str
        Gives title of the plot.
 
    Returns
    -------
    symbols : 1D array of strs
        An array containing strings of either an N or V, corresponding to normal or abnormal event.
    trial_time : 1D array of floats
        Gives a 1D array of 250 times between -0.5 and 0.5 seconds.
    mean_trial_signal : 2D array of size 2X250 that gives mean signal values at each sample.
        Gives a 2D array of size 2x250 giving the means of each column for the n and v trials.
 
    """
    
    trial_time = np.arange(-(trial_duration_seconds/2),trial_duration_seconds/2,1/fs)
    symbols = np.unique(label_symbols)
    
    # Setting plot
    plt.figure(3, figsize=(10, 6), clear=True)
    plt.xlabel('Time (sec)')
    plt.ylabel(f'ECG Voltage ({units})')
    plt.title(title)
    plt.grid(True)
    
    # Initialize mean trial signal 2D array
    mean_trial_signal = np.zeros((2,fs), dtype=float)
    
    # Count corresponds to row in mean_trial_signal
    mean_count = 0
    
    # Loop through N and V label
    for label in symbols:
        event_samples_label = label_samples[label_symbols == label]
        start_indices = np.zeros(len(event_samples_label), dtype=int)
        
        for i, event_index in enumerate(event_samples_label):
            start_indices[i] = max(event_index - fs/2, 0) #fs/2 = 125
            
        trials = extract_trials(signal_voltage, start_indices, fs) #fs = 250
        
        # Get means and std deviations for each label in symbols
        means = np.mean(trials, axis=0)
        stds = np.std(trials, axis=0)
        mean_trial_signal[mean_count,:] = means
        mean_count += 1
        
        # Plot trials mean and shaded region for standard deviation
        plt.plot(trial_time, means, label=f'{label} Trials Mean',linewidth=1.5)
        plt.fill_between(trial_time, means - stds, means + stds, alpha=0.2, label=f'{label} Trials ±1 Std Dev')
        
    plt.legend()
    
    return symbols, trial_time, mean_trial_signal

#%% Part 6: Save arrays and plots

def save_means(symbols, trial_time, mean_trial_signal, out_filename):
    """
    Parameters
    ----------
    symbols : 1D array of strings of size 2x1
        Gives either an N or V.
    trial_time : 1D array of floats, size 250x1
        Gives 250 time samples between -0.5 and 0.5 seconds.
    mean_trial_signal : 2D array of floats of size 2X250.
        Gives the means of each column for the n and v trials.
    out_filename : string
        Gives the name of the output file created.
 
    Returns
    -------
    None.
 
    """
    
    # Save data as npz file
    np.savez(out_filename, symbols=symbols, trial_time=trial_time, mean_trial_signal=mean_trial_signal)
    
    
    
