import matplotlib.pyplot as plt 
import numpy as np 
import os 
import plotly.express as px 
import argparse 


singulars = [20, 40, 80, 160, 320, 640, 1280, 2560, 3355]

def plot_images(force):
    data = np.loadtxt('dog_bw_data.dat')
    plt.gray()
    plt.imshow(data)
    plt.savefig('dog_original.png')

    print('Generating compressed images')
    for i, num in enumerate(singulars):
        if not os.path.isfile(f'compressed_{num}.png') or force:
            print(f'Reading in {i}')
            data = np.loadtxt(f'compr {i+1}')
            plt.gray()
            plt.imshow(data)
            plt.savefig(f'compressed_{num}.png')
        else:
            print(f'Reconstruction with {num = }th plot exists, continuing...')

def calculate_errors(force):
    data = np.loadtxt('dog_bw_data.dat')
    
    print('Generating error plot')
    errors = np.loadtxt('frobenius_errors.dat')
    fig = px.line(x=singulars, y=errors)
    fig.write_image('errors.png', scale=2)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--force',
        required=False,
        action='store_true',
        default=False
    )

    force = parser.parse_args().force 
