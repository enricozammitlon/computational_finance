import matplotlib.pyplot as plt
import csv
import numpy as np


def firstTask():
    analyticData = {'S': [], 'analytic_pi_0': [], 'analytic_pi_T': []}
    with open('outputs/analytic.task.2.1.csv', newline='\n') as csvfile:
        reader = csv.DictReader(csvfile, fieldnames=[
            'S', 'analytic_pi_0', 'analytic_pi_T'], quoting=csv.QUOTE_NONNUMERIC)
        for row in reader:
            analyticData['S'].append(row['S'])
            analyticData['analytic_pi_0'].append(row['analytic_pi_0'])
            analyticData['analytic_pi_T'].append(row['analytic_pi_T'])
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        plt.grid(True)
        ax.plot(analyticData['S'], analyticData['analytic_pi_0'],
                label=r't=0', color='red', linewidth=2)
        ax.plot(analyticData['S'], analyticData['analytic_pi_T'],
                label=r't=T', color='blue', linewidth=2)
        plt.xlabel('S')
        plt.ylabel('V(S,t)')
        plt.legend()
        plt.savefig('Solution/analytic_values.png',
                    bbox_inches='tight', pad_inches=0.25)

    files = ['normal.task.2.1', 'antithetic.task.2.1', 'momentmatch.task.2.1']
    allData = {'70': {'normal': {}, 'antithetic': {}, 'momentmatch': {}},
               '100': {'normal': {}, 'antithetic': {}, 'momentmatch': {}}}
    for file in files:
        counter = 0
        with open('outputs/'+file+'.csv', newline='\n') as csvfile:
            # N, S0, Montecarlo Avrg PI, Lower confidence limit,Upper confidence limit, Analytic PI
            reader = csv.DictReader(csvfile, fieldnames=[
                                    'time_taken', 'n', 's0', 'montecarlo_pi', 'std', 'lcl', 'ucl', 'analytic_pi'], quoting=csv.QUOTE_NONNUMERIC)
            currentFileData = {'time_taken': [], 'n': [], 's0': [], 'montecarlo_pi': [
            ], 'std': [], 'error_bars': [[], []], 'analytic_pi': []}
            for row in reader:
                counter += 1
                if(counter == 101):
                    allData[str(int(currentFileData['s0'][0]))
                            ][file.split(sep='.')[0]] = currentFileData
                    currentFileData = {'time_taken': [], 'n': [], 's0': [], 'montecarlo_pi': [
                    ], 'std': [], 'error_bars': [[], []], 'analytic_pi': []}
                currentFileData['n'].append(row['n']/1000)
                currentFileData['time_taken'].append(row['time_taken'])
                currentFileData['s0'].append(row['s0'])
                currentFileData['montecarlo_pi'].append(row['montecarlo_pi'])
                currentFileData['std'].append(
                    100*row['std']/row['montecarlo_pi'])
                currentFileData['error_bars'][0].append(
                    row['montecarlo_pi']-2*row['std'])
                currentFileData['error_bars'][1].append(
                    row['montecarlo_pi']+2*row['std'])
                currentFileData['analytic_pi'].append(row['analytic_pi'])
            allData[str(int(currentFileData['s0'][0]))
                    ][file.split(sep='.')[0]] = currentFileData

    locations = {'70_': '', }
    for dataKey, dataFrame in allData.items():
        for subKey, subFrame in dataFrame.items():
            fig, ax1 = plt.subplots()
            plt.grid(True)
            ax2 = ax1.twinx()
            ax2.tick_params('both', labelsize=14)
            ax1.tick_params('both', labelsize=14)

            ax1.set_xlabel('N / Thousands', fontsize=16)
            ax1.set_ylabel('V($S_0$='+dataKey+',t=0)', fontsize=16)
            ax2.set_ylabel('Standard Error (%)', fontsize=16)
            # p3 = np.poly1d(np.polyfit(subFrame['n'], subFrame['std'], 5))
            ax2.plot(subFrame['n'], subFrame['std'], label=r'(Standard Error)',
                     ls='-', color='red', linewidth=2)
            # This is about confidence intervals
            n = np.array(subFrame['n'], dtype=np.float64)
            mean = np.array(subFrame['montecarlo_pi'], dtype=np.float64)
            lower = np.array(subFrame['error_bars'][0], dtype=np.float64)
            upper = np.array(subFrame['error_bars'][1], dtype=np.float64)
            analytic = np.array(subFrame['analytic_pi'], dtype=np.float64)
            ax1.plot(n, mean,
                     label=r'Montecarlo Average', color='blue', alpha=0.7, linewidth=1)
            ax1.plot(n, analytic,
                     label=r'Analytic Value', ls='dashed', color='green', linewidth=3)
            ax1.set_xlim(1, subFrame['n'][-1])
            handles, labels = [
                (a + b) for a, b in zip(ax1.get_legend_handles_labels(), ax2.get_legend_handles_labels())]
            plt.savefig('Solution/confidence_'+dataKey+'_'+subKey+'.png',
                        bbox_inches='tight', pad_inches=0.25)

            axe = plt.axes(frameon=False)
            axe.figure.set_size_inches(8, 1)
            axe.legend(handles, labels, ncol=3, loc='center',
                       mode="expand", fancybox=False, framealpha=0.0)
            axe.xaxis.set_visible(False)
            axe.yaxis.set_visible(False)
            plt.savefig('Solution/legend.png',
                        bbox_inches='tight', pad_inches=0)

    labels = {'antithetic': 'Antithetic Variables',
              'normal': 'Basic Montecarlo', 'momentmatch': 'Moment Matching'}
    colors = {'antithetic': 'red',
              'normal': 'blue', 'momentmatch': 'green'}

    fig_70 = plt.figure()
    fig_100 = plt.figure()
    ax_70 = fig_70.add_subplot(1, 1, 1)
    ax_100 = fig_100.add_subplot(1, 1, 1)

    for dataKey, dataFrame in allData.items():
        for subDataKey, subDataFrame in dataFrame.items():
            # This is about timing efficiency
            ax_70.plot(allData["70"][subDataKey]['n'], allData["70"][subDataKey]['time_taken'],
                       label=labels[subDataKey], color=colors[subDataKey], linewidth=2)

            ax_100.plot(allData["100"][subDataKey]['n'], allData["100"][subDataKey]['time_taken'],
                        label=labels[subDataKey], color=colors[subDataKey], linewidth=2)
        break

    plt.legend()
    plt.grid(True)
    plt.xlabel('N')
    plt.ylabel('Time / ms')
    fig_70.savefig('Solution/timing_efficiency_70.png',
                   bbox_inches='tight', pad_inches=0.25)
    fig_100.savefig('Solution/timing_efficiency_100.png',
                    bbox_inches='tight', pad_inches=0.25)

# Fixed N,K and vary M


def vary_M():
    files = ['paths.task.2.2']
    allData = {'K': []}
    value_of_M = 0

    for file in files:
        counter = 0
        with open('outputs/'+file+'.csv', newline='\n') as csvfile:
            # N, S0, Montecarlo Avrg PI, Lower confidence limit,Upper confidence limit, Analytic PI
            reader = csv.DictReader(csvfile, fieldnames=[
                                    'time_taken', 'n', 'k', 'm', 'montecarlo_pi', 'std', 'lcl', 'ucl'], quoting=csv.QUOTE_NONNUMERIC)
            currentFileData = {'time_taken': [], 'n': [], 'k': [], 'm': [], 'montecarlo_pi': [
            ], 'std': [], 'error_bars': [[], []]}
            found = False
            for row in reader:
                counter += 1
                if(counter == 46):
                    counter = 1
                    allData['K'].append(currentFileData)
                    currentFileData = {'time_taken': [], 'n': [], 'k': [], 'm': [], 'montecarlo_pi': [
                    ], 'std': [], 'error_bars': [[], []]}
                currentFileData['n'].append(row['n']/1000)
                currentFileData['time_taken'].append(row['time_taken'])
                currentFileData['k'].append(row['k'])
                currentFileData['montecarlo_pi'].append(row['montecarlo_pi'])
                currentFileData['std'].append(
                    100*row['std']/row['montecarlo_pi'])
                currentFileData['error_bars'][0].append(
                    row['montecarlo_pi']-2*row['std'])
                currentFileData['error_bars'][1].append(
                    row['montecarlo_pi']+2*row['std'])
                currentFileData['m'].append(row['m'])
                if(found != True and row['time_taken'] >= 10000):
                    found = True
                    value_of_M = row['m']
            allData['K'].append(currentFileData)

    fig, ax1 = plt.subplots()
    plt.grid(True)
    ax1.tick_params('both', labelsize=14)

    ax1.set_xlabel('Paths', fontsize=16)
    ax1.set_ylabel('Option Price', fontsize=16)
    colors = {'20': 'blue', '50': 'green', '100': 'red', '150': 'purple'}
    for dataKey, dataFrame in allData.items():
        for subFrame in dataFrame:
            # This is about number of paths
            m = np.array(subFrame['m'], dtype=np.float64)
            mean = np.array(subFrame['montecarlo_pi'], dtype=np.float64)
            ax1.plot(m, mean,
                     label=r'$K$='+str(subFrame['k'][0]), color=colors[str(int(subFrame['k'][0]))], alpha=1, linewidth=1)
        handles, labels = ax1.get_legend_handles_labels()
        fig.legend(handles, labels)
        fig.savefig('Solution/task_2_2_plot_1.png',
                    bbox_inches='tight', pad_inches=0.25)
    fig, ax2 = plt.subplots()
    plt.grid(True)
    ax2.tick_params('both', labelsize=14)

    ax2.set_xlabel('Paths', fontsize=16)
    ax2.set_ylabel('Standard Error (%)', fontsize=16)
    ax2.axvline(x=value_of_M, ymin=0, ls='--',
                color='purple', label=r"10s Computation")
    found = False
    for dataKey, dataFrame in allData.items():
        for subFrame in dataFrame:
            # This is about standard error
            ax2.plot(subFrame['m'], subFrame['std'],
                     label=r'$K$='+str(subFrame['k'][0]), color=colors[str(int(subFrame['k'][0]))], linewidth=1)

        handles, labels = ax2.get_legend_handles_labels()
        fig.legend(handles, labels)
        fig.savefig('Solution/task_2_2_plot_2.png',
                    bbox_inches='tight', pad_inches=0.25)

# Fixed N,M and vary K


def vary_K():
    files = ['k.task.2.2']
    allData = {}
    for file in files:
        counter = 0
        with open('outputs/'+file+'.csv', newline='\n') as csvfile:
            # N, S0, Montecarlo Avrg PI, Lower confidence limit,Upper confidence limit, Analytic PI
            reader = csv.DictReader(csvfile, fieldnames=[
                                    'time_taken', 'n', 'k', 'm', 'montecarlo_pi', 'std', 'lcl', 'ucl'], quoting=csv.QUOTE_NONNUMERIC)
            allData = {'time_taken': [], 'n': [], 'k': [], 'm': [], 'montecarlo_pi': [
            ], 'std': [], 'error_bars': [[], []]}
            for row in reader:
                allData['n'].append(row['n']/1000)
                allData['time_taken'].append(row['time_taken'])
                allData['k'].append(row['k'])
                allData['montecarlo_pi'].append(row['montecarlo_pi'])
                allData['std'].append(
                    100*row['std']/row['montecarlo_pi'])
                allData['error_bars'][0].append(
                    row['montecarlo_pi']-2*row['std'])
                allData['error_bars'][1].append(
                    row['montecarlo_pi']+2*row['std'])
                allData['m'].append(row['m'])

    fig, ax1 = plt.subplots()
    plt.grid(True)
    ax1.tick_params('both', labelsize=14)

    ax1.set_xlabel('K', fontsize=16)
    ax1.set_ylabel('Option Price', fontsize=16)

    # This is about number of paths
    k = np.array(allData['k'], dtype=np.float64)
    mean = np.array(allData['montecarlo_pi'], dtype=np.float64)
    ax1.scatter(k, mean)
    #ax1.set_xlim(20, allData['k'][-1])
    fig.savefig('Solution/task_2_2_plot_3.png',
                bbox_inches='tight', pad_inches=0.25)

    fig, ax2 = plt.subplots()
    plt.grid(True)
    ax2.tick_params('both', labelsize=14)
    ax2.set_xlabel('K', fontsize=16)
    ax2.set_ylabel('Standard Error (%)', fontsize=16)
    #ax2.set_xlim(20, allData['k'][-1])
    # This is about standard error
    ax2.plot(allData['k'], allData['std'], linewidth=1)

    handles, labels = ax2.get_legend_handles_labels()
    fig.legend(handles, labels)
    fig.savefig('Solution/task_2_2_plot_4.png',
                bbox_inches='tight', pad_inches=0.25)


def vary_N():
    files = ['n.task.2.2']
    allData = {'K': []}
    value_of_N = 0
    for file in files:
        counter = 0
        found = False
        with open('outputs/'+file+'.csv', newline='\n') as csvfile:
            # N, S0, Montecarlo Avrg PI, Lower confidence limit,Upper confidence limit, Analytic PI
            reader = csv.DictReader(csvfile, fieldnames=[
                                    'time_taken', 'n', 'k', 'm', 'montecarlo_pi', 'std', 'lcl', 'ucl'], quoting=csv.QUOTE_NONNUMERIC)
            currentFileData = {'time_taken': [], 'n': [], 'k': [], 'm': [], 'montecarlo_pi': [
            ], 'std': [], 'error_bars': [[], []]}
            for row in reader:
                counter += 1
                if(counter == 51 or counter == 92):
                    allData['K'].append(currentFileData)
                    currentFileData = {'time_taken': [], 'n': [], 'k': [], 'm': [], 'montecarlo_pi': [
                    ], 'std': [], 'error_bars': [[], []]}
                currentFileData['n'].append(row['n']/1000)
                currentFileData['time_taken'].append(row['time_taken'])
                currentFileData['k'].append(row['k'])
                currentFileData['montecarlo_pi'].append(row['montecarlo_pi'])
                currentFileData['std'].append(
                    100*row['std']/row['montecarlo_pi'])
                currentFileData['error_bars'][0].append(
                    row['montecarlo_pi']-2*row['std'])
                currentFileData['error_bars'][1].append(
                    row['montecarlo_pi']+2*row['std'])
                currentFileData['m'].append(row['m'])
                if(found != True and row['time_taken'] >= 10000):
                    found = True
                    value_of_N = row['n']/1000
            allData['K'].append(currentFileData)

    fig, ax1 = plt.subplots()
    plt.grid(True)
    ax1.tick_params('both', labelsize=14)

    ax1.set_xlabel('N/Thousands', fontsize=16)
    ax1.set_ylabel('Option Price', fontsize=16)
    colors = {'20': 'blue', '50': 'green', '70': 'red'}
    for dataKey, dataFrame in allData.items():
        for subFrame in dataFrame:
            # This is about number of paths
            n = np.array(subFrame['n'], dtype=np.float64)
            mean = np.array(subFrame['montecarlo_pi'], dtype=np.float64)
            ax1.plot(n, mean,
                     label=r'$K$='+str(subFrame['k'][0]), color=colors[str(int(subFrame['k'][0]))], alpha=1, linewidth=1)
        handles, labels = ax1.get_legend_handles_labels()
        fig.legend(handles, labels)
        fig.savefig('Solution/task_2_2_plot_5.png',
                    bbox_inches='tight', pad_inches=0.25)
    fig, ax2 = plt.subplots()
    plt.grid(True)
    ax2.tick_params('both', labelsize=14)

    ax2.set_xlabel('N/Thousands', fontsize=16)
    ax2.set_ylabel('Standard Error (%)', fontsize=16)
    ax2.axvline(x=value_of_N, ymin=0, ls='--',
                color='purple', label=r"10s Computation")
    for dataKey, dataFrame in allData.items():
        for subFrame in dataFrame:
            # This is about standard error
            ax2.plot(subFrame['n'], subFrame['std'],
                     label=r'$K$='+str(subFrame['k'][0]), color=colors[str(int(subFrame['k'][0]))], linewidth=1)

        handles, labels = ax2.get_legend_handles_labels()
        fig.legend(handles, labels)
        fig.savefig('Solution/task_2_2_plot_6.png',
                    bbox_inches='tight', pad_inches=0.25)


vary_M()
vary_K()
vary_N()
