import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.serif'] = ['Linux Biolinum', 'Arial']
mpl.rcParams['font.sans-serif'] = ['Linux Libertine', 'Times New Roman']
mpl.rcParams['font.size'] = 14
mpl.rcParams['mathtext.fontset'] = 'custom'

SCHEMES = ["drape", "twisting", "locking", "ball"]#, "armadillo", "needles", "umbrella", "can (no mandrel)", "can (with mandrel)"]
COLORS = ['tab:blue', 'tab:green', 'tab:orange']
NAMES = {"hang":[["hang10x10","hang20x20","hang30x30"],["new4k","new16k","new64k"],["phang10","phang20","phang30"]],
        "wrinkle":[["wrinkle10x10","wrinkle15x15","wrinkle20x20"],["wrinkle32k","wrinkle65k","wrinkle130k"],["pseudo10","pseudo15","pseudo20"]]}
OUTPUTS = {"hang":[["BHEM10$\\times$10","BHEM20$\\times$20","BHEM30$\\times$30"],["FEM 4k","FEM 16k","FEM 64k"],["pseudo10","pseudo20","pseudo30"]],
        "wrinkle":[["BHEM10$\\times$10","BHEM15$\\times$15","BHEM20$\\times$20"],["FEM 32k","FEM 65k","FEM 130k"],["pseudo10","pseudo15","pseudo20"]]}
XLIM = {"hang":60, "wrinkle":50}

def plot_iter(filename):
    fig, ax = plt.subplots()
    for i, name in enumerate(NAMES[filename][0]):
        data = np.loadtxt(name+".txt", dtype=np.float32,delimiter='  ') 
        # print(data)
        errors = data[:,0]
        times = data[:,1]
        ax.semilogy(range(len(errors)), np.array(errors), label=OUTPUTS[filename][0][i], color=COLORS[i])
    for i, name in enumerate(NAMES[filename][1]):
        data = np.loadtxt(name+".txt", dtype=np.float32,delimiter='  ') 
        # print(data)
        errors = data[:,0]
        times = data[:,1]
        ax.semilogy(range(len(errors)), np.array(errors), label=OUTPUTS[filename][1][i], color=COLORS[i], linestyle='dashed')
    for i, name in enumerate(NAMES[filename][2]):
        data = np.loadtxt(name+".txt", dtype=np.float32,delimiter='  ') 
        # print(data)
        errors = data[:,0]
        times = data[:,1]
        ax.semilogy(range(len(errors)), np.array(errors), label=OUTPUTS[filename][2][i], color=COLORS[i], linestyle='dotted')
    ax.legend(prop={'family': 'serif', 'size': 10}, loc='upper right')
    ax.set_xlabel('Number of Iterations', fontfamily='sans-serif')
    # ax.set_xlabel('$\Delta x\,$/m', fontfamily='sans-serif')
    ax.set_ylabel('Relative Residual', fontfamily='serif')
    # ax.set_xticks([0.01, 0.02, 0.05, 0.1])
    plt.ylim(1e-10,1e0)
    plt.xlim(0,XLIM[filename])
    ax.grid(color='tab:grey', alpha=0.5, linestyle='dashed', linewidth=0.5)

    fig.savefig(filename + '-iter.eps', dpi=300, bbox_inches='tight', pad_inches=0.0)
    fig.savefig(filename + '-iter.png', dpi=300)

def plot_time(filename):
    fig, ax = plt.subplots()
    for i, name in enumerate(NAMES[filename][0]):
        data = np.loadtxt(name+".txt", dtype=np.float32,delimiter='  ') 
        # print(data)
        errors = data[:,0]
        times = np.append([0],data[:-1,1])
        ax.semilogy(np.array(times), np.array(errors), label=OUTPUTS[filename][0][i], color=COLORS[i])
    for i, name in enumerate(NAMES[filename][1]):
        data = np.loadtxt(name+".txt", dtype=np.float32,delimiter='  ') 
        # print(data)
        errors = data[:,0]
        times = np.append([0],data[:-1,1])
        ax.semilogy(np.array(times), np.array(errors), label=OUTPUTS[filename][1][i], color=COLORS[i], linestyle='dashed')
    for i, name in enumerate(NAMES[filename][2]):
        data = np.loadtxt(name+".txt", dtype=np.float32,delimiter='  ') 
        # print(data)
        errors = data[:,0]
        times = np.append([0],data[:-1,1])
        ax.semilogy(np.array(times), np.array(errors), label=OUTPUTS[filename][2][i], color=COLORS[i], linestyle='dotted')
    ax.legend(prop={'family': 'serif', 'size': 10}, loc='upper right')
    ax.set_xlabel('Time Cost (s)', fontfamily='sans-serif')
    # ax.set_xlabel('$\Delta x\,$/m', fontfamily='sans-serif')
    ax.set_ylabel('Relative Residual', fontfamily='serif')
    # ax.set_xticks([0.01, 0.02, 0.05, 0.1])
    plt.ylim(1e-10,1e0)
    plt.xlim(0,100)
    ax.grid(color='tab:grey', alpha=0.5, linestyle='dashed', linewidth=0.5)

    fig.savefig(filename + '-time.eps', dpi=300, bbox_inches='tight', pad_inches=0.0)
    fig.savefig(filename + '-time.png', dpi=300)


if __name__ == '__main__':
    plot_iter("hang")
    plot_time("hang")