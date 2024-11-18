from scipy.ndimage import map_coordinates
import matplotlib.animation as ani
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import pandas as pd
import numpy as np
import os

# physics constants
c0 = 299792458  # m/s
G0 = 6.67430e-11  # m^3/kg/s^2
Msun = 0.989e30  # kg
mpc = 2.08567758e22  # m
# NGC 1300 Info
BG_DISTANCE = 16.5  # Mpc
# file paths
DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
(_, _, DATA_FILES) = next(os.walk(DATA_DIR))

def read(file: str) -> tuple[2]:
    with fits.open(file) as h:
        header, data = h[0].header, h[0].data
    return (header, data)

class BH:
    def __init__(self, mass: float, distance: float, x: float, y: float):
        self.mass = mass            # Solar masses
        self.distance = distance    # Mpc
        self.x = x                  # x-coordinate(pixel)
        self.y = y                  # y-coordinate(pixel)
    
    def __str__(self) -> str:
        return f'BH(mass={self.mass}, distance={self.distance}, x={self.x}, y={self.y})'

class Plotter:
    def __init__(self, data: tuple):
        self.header, self.data = data
        xpixelsz = self.metadata.loc[self.metadata['Key'] == 'XPIXELSZ'].iloc[0]['Value']
        ypixelsz = self.metadata.loc[self.metadata['Key'] == 'YPIXELSZ'].iloc[0]['Value']
        pltscale = self.metadata.loc[self.metadata['Key'] == 'PLTSCALE'].iloc[0]['Value']
        # print(f'XPIXELSZ: {xpixelsz} um/px, YPIXELSZ: {ypixelsz} um/px, PLTSCALE: {pltscale} arcsec/mm')
        self.xscale = xpixelsz * pltscale / 1000 * np.pi / 180 / 3600   # rad/px
        self.yscale = ypixelsz * pltscale / 1000 * np.pi / 180 / 3600   # rad/px

    @property
    def metadata(self):
        return pd.DataFrame(self.header.items(), columns=['Key', 'Value'])
    
    def plot(self, output: str = None, show: bool = True):
        fig, ax = plt.subplots(figsize=(12, 12))
        ax.imshow(self.data, cmap='gray')
        if show:
            plt.show()
        if output:
            fig.savefig(output)
    
    def apply_lensing(self, bh: BH):
        # Constants
        D_l = bh.distance * mpc          # Lens distance in meters
        D_s = BG_DISTANCE * mpc          # Source distance in meters
        D_ls = D_s - D_l                 # Distance between lens and source
        M = bh.mass * Msun               # Mass in kg

        # Pixel grid
        y_indices, x_indices = np.indices(self.data.shape)
        x = (x_indices - bh.x) * self.xscale  # x in radians
        y = (y_indices - bh.y) * self.yscale  # y in radians
        theta = np.sqrt(x**2 + y**2)          # Angular position
        epsilon = 1e-6
        theta = np.where(theta == 0, epsilon, theta)    # replace zero with eps

        # Mapping
        alpha = (4 * G0 * M) / (c0**2 * theta * D_l)    # Deflection angle
        theta_s = theta - alpha * (D_ls / D_s)          # Source position
        theta_s_over_theta = theta_s / theta            # Magnification
        x_s = x * theta_s_over_theta                    # x in radians
        y_s = y * theta_s_over_theta                    # y in radians
        x_s_px = x_s / self.xscale + bh.x               # convert back to x_px
        y_s_px = y_s / self.yscale + bh.y               # convert back to y_px

        # Directly output
        lensed = np.zeros(self.data.shape)
        valid = (x_s_px >= 0) & (x_s_px < self.data.shape[1]) & (y_s_px >= 0) & (y_s_px < self.data.shape[0])
        x_s_px = x_s_px[valid]
        y_s_px = y_s_px[valid]
        lensed[y_indices[valid], x_indices[valid]] = self.data[y_s_px.astype(int), x_s_px.astype(int)]

        # Interpolate the image
        # coords = np.array([y_s_px.flatten(), x_s_px.flatten()])
        # lensed_flat = map_coordinates(self.data, coords, order=1, mode='nearest')
        # lensed = lensed_flat.reshape(self.data.shape)
        return lensed

def animation(units: list, images: list, itval: int, output: str = None):
    fig, ax = plt.subplots(figsize=(12, 12))
    # ax.axis('off')
    im = ax.imshow(images[0], cmap='gray')
    def update(i):
        im.set_data(images[i])
        ax.set_title(f'variable: {units[i]:.2e} units')
        return im,
    anime = ani.FuncAnimation(fig, update, frames=len(images), interval=itval)
    if output:
        anime.save(os.path.join(DATA_DIR, output + ".gif"), writer='imagemagick')
    
if __name__ == "__main__":
    for file in DATA_FILES:
        if not file.endswith('.fits'):
            continue
        newPlotter = Plotter(read(os.path.join(DATA_DIR, file)))
        newPlotter.plot(output=os.path.join(DATA_DIR, file.replace('.fits', '.png')), show=False)

        # 1. change mass
        masses = np.logspace(11, 14, 100)
        images = [newPlotter.apply_lensing(BH(mass, 10, 385, 510)) for mass in masses]
        animation(masses, images, 100, output=file.replace('.fits', '_lensing1'))

        # 2. change distance
        dists = np.linspace(1, 16, 100)
        images = [newPlotter.apply_lensing(BH(1e13, dist, 300, 300)) for dist in dists]
        animation(dists, images, 100, output=file.replace('.fits', '_lensing2'))

        # 3. change position
        xs = np.linspace(100, 500, 101, dtype=int)
        images = [newPlotter.apply_lensing(BH(1e13, 10, x, 300)) for x in xs]
        animation(xs, images, 100, output=file.replace('.fits', '_lensing3'))