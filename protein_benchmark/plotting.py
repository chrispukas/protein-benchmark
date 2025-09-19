import numpy as np
import matplotlib.pyplot as plt




"""
PLOTTING STACKED BAR CHART UTILITY
"""

def plot_stacked_bar_chart(categorised_data_dict, 
                           colours=['#e377c2', '#7bac4e', '#d4c04c', '#4c70c4', '#8c564b', '#e7c377', '#c7c7c7'], 
                           title="", 
                           legend_title="DockQ Categories",
                           y_label="Percentage of Models (%)", 
                           remove_low_cutoff=True,
                           fixed_max_y=None,
                           figsize: tuple = (4,6)):
    plt.figure(figsize=figsize)

    group_labels = list(categorised_data_dict.keys())
    categorised_data_values = list(categorised_data_dict.values())
    sub_categories = list(categorised_data_values[0].keys())

    for i, data in enumerate(categorised_data_values):
        total = sum(len(data[cat]) for cat in data) if fixed_max_y is None else fixed_max_y
        bottom = 0
        for j, cat in enumerate(reversed(sub_categories)):
            if remove_low_cutoff and cat == "DockQ<0.23":
                continue
            value = (len(data.get(cat, [])) / total) * 100 if total > 0 else 0
            plt.bar(group_labels[i],
                    value,
                    bottom=bottom,
                    color=colours[j % len(colours)],
                    label=cat if i == 0 else "_nolegend_")
            bottom += value

    plt.yticks(np.arange(0, 101, 10))
    plt.ylabel(y_label)
    plt.ylim(0, 100)
    
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), title=legend_title)
    plt.title(title)
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    plt.show()


"""
PLOTTING SIDED BAR CHART UTILITY
"""

def plot_sided_bar_chart(categorised_data_dict, colours=['#4c70c4', '#7bac4e', '#d4c04c', '#e377c2', '#8c564b', '#e7c377', '#c7c7c7'], title="", legend_title="DockQ Categories", remove_low_cutoff=True):
    plt.figure(figsize=(4, 6))

    group_labels = list(categorised_data_dict.keys())
    values = list(categorised_data_dict.values())

    for i, data in enumerate(values):
        plt.bar(group_labels[i],
                data,
                0.8,
                color=colours[i % len(colours)],)




    plt.yticks(np.arange(0.0, 1.1, 0.1))
    plt.ylabel('DockQ Score')
    plt.ylim(0, 1.0)
    
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), title=legend_title)
    plt.title(title)
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    plt.show()



"""
PLOTTING VIOLIN PLOT UTILITY
"""

def plot_violin_plot(categorized_data_dict, colours=['#4c70c4', '#7bac4e', '#d4c04c', '#e377c2', '#8c564b', '#e7c377', '#c7c7c7'], title=""):
    group_labels = list(categorized_data_dict.keys())
    categorized_data = list(categorized_data_dict.values())

    plt.figure(figsize=(4, 6))

    print(categorized_data_dict)

    for i, data in enumerate(categorized_data):
        print(data)
        plt.violinplot(data, positions=[i], showmeans=True, showmedians=True, widths=0.8, 
                        points=100, bw_method='scott', vert=True)

    plt.yticks(np.arange(0, 1.1, 0.1))
    plt.ylim(0, 1)
    plt.xticks(np.arange(len(group_labels)), group_labels)
    plt.ylabel('DockQ Scores')
    plt.title(title)
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    plt.show()


"""
PLOTTING SCATTER PLOT UTILITY
"""

def plot_scatter_plot(x_data, y_data, x_label="X-axis", y_label="Y-axis", title="Scatter Plot", scatter_colour='blue', line_colour='black', marker='o', xlim=None, ylim=None):
    plt.figure(figsize=(6, 4))
    plt.scatter(x_data, y_data, color=scatter_colour, marker=marker)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)

    ax = plt.gca()
    ax.axline(xy1=(0, 0), slope=1, color=line_colour, linestyle='--', linewidth=1.5)

    plt.title(title)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.show()