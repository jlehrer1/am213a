import matplotlib.pyplot as plt 
import numpy as np 
import os 
import plotly.express as px 
import plotly.graph_objects as go 
import argparse

def plot_images(singulars, force):
    data = np.loadtxt('dog_bw_data.dat')
    plt.gray()
    plt.imshow(data)
    plt.savefig('dog_original.png')

    print('Generating compressed images')
    for i, num in enumerate(singulars):
        if not os.path.isfile(f'images/compressed_{num}.png') or force:
            print(f'Reading in {i}')
            data = np.loadtxt(f'images/compr {i+1}')
            plt.gray()
            plt.axis('off')
            plt.imshow(data)
            plt.savefig(f'compressed_{num}.png')
        else:
            print(f'Reconstruction with {num = }th plot exists, continuing...')

def plot_parta_errors(singulars, force):
    filename = 'errors.png'

    if not os.path.isfile(filename) or force:
        errors = np.loadtxt('frobenius_errors.dat')
        fig = px.line(x=singulars, y=errors)
        fig.update_layout(
            xaxis_title='Number of Singular Values',
            yaxis_title='Average Frobenius Errors',
            title='Frobenius Compression Error vs Number of Singular Values'
        )

        fig.write_image(filename, scale=3)
    else:
        print('Error plots exist, continuing...')

def plot_individual(errfiles, plotname, force, title):
    ds = [2, 5, 10, 100, 1000]

    if not os.path.isfile(plotname) or force:
        data = []
        for errfile, d in zip(errfiles, ds):
            errs = np.loadtxt(errfile)[0:20]

            data.append((
                go.Scatter(
                    x=list(range(len(errs))),
                    y=np.log(errs),
                    mode='lines',
                    name=f'{d =}'
                )
            ))

        fig = go.Figure(
            data=data,
            layout=go.Layout(
                xaxis=dict(title='# Iterations'),
                yaxis=dict(title='# Norm Error of 'r'$Ax^{(k)} - b$'),
                title=title,
            )
        )

        fig.write_image(plotname, scale=2)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--force',
        required=False,
        action='store_true',
        default=False,
        help='Re-calculate plots, even if they exist'
    )

    parser.add_argument(
        '--parta',
        required=False,
        action='store_true',
        default=False,
        help='If passed, only generate plots for part a. Otherwise, generate everything.'
    )

    parser.add_argument(
        '--partb',
        required=False,
        action='store_true',
        default=False,
        help='If passed, only generate plots for part b. Otherwise, generate everything.'
    )

    args = parser.parse_args()
    force, parta, partb = args.force, args.parta, args.partb 

    singulars = [20, 40, 80, 160, 320, 640, 1280, 2560, 3355]

    gs_errors = [f'errors/GS_er {i}' for i in range(1, 6)]
    gj_errors = [f'errors/GJ_er {i}' for i in range(1, 6)]

    gs_errfile = 'images/gauss_seidel_error.png'
    gj_errfile = 'images/gauss_jacobi_error.png'

    if parta:
        plot_images(singulars, force)
        plot_parta_errors(singulars, force)
        exit()
    if partb:
        plot_individual(gs_errors, gs_errfile, force, 'Gauss-Seidel, Log 2-Norm Error')
        plot_individual(gj_errors, gj_errfile, force, 'Gauss-Jacobi, Log 2-Norm Error')
        exit()

    plot_images(singulars, force)
    plot_parta_errors(singulars, force)
    plot_individual(gs_errors, gs_errfile, force, 'Gauss-Seidel, Log 2-Norm Error')
    plot_individual(gj_errors, gj_errfile, force, 'Gauss-Jacobi, Log 2-Norm Error')
