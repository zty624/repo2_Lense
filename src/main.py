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
        rs = 2 * G0 * bh.mass * Msun / c0**2    # Schwarzschild radius
        # center the black hole at (0, 0)
        ypx0, xpx0 = np.indices(self.data.shape)
        xpx = xpx0 - bh.x
        ypx = ypx0 - bh.y

        theta = np.sqrt((xpx * self.xscale)**2 + (ypx * self.yscale)**2)
        b = theta * bh.distance * mpc + 1e5
        alpha = 2 * rs/ (b)
        r = np.sqrt(xpx ** 2 + ypx ** 2) + 1e5

        new_xpx = xpx - alpha / r * xpx * (BG_DISTANCE - bh.distance) / BG_DISTANCE / self.xscale
        new_ypx = ypx - alpha / r * ypx * (BG_DISTANCE - bh.distance) / BG_DISTANCE / self.yscale
        new_xpx = np.clip(new_xpx.astype(int) + bh.x, 0, self.data.shape[1] - 1)
        new_ypx = np.clip(new_ypx.astype(int) + bh.y, 0, self.data.shape[0] - 1)

        lensed = np.zeros_like(self.data)
        cut = (new_xpx >= 0) & (new_xpx < self.data.shape[1]) & (new_ypx >= 0) & (new_ypx < self.data.shape[0])
        lensed[new_ypx[cut], new_xpx[cut]] = self.data[ypx0[cut], xpx0[cut]]
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

        # # 1. change mass
        # masses = np.logspace(15, 16.5, 100)
        # images = [newPlotter.apply_lensing(BH(mass, 10, 522, 424)) for mass in masses]
        masses = np.logspace(15, 16.5, 100)
        images = [newPlotter.apply_lensing(BH(mass, 10, 385, 510)) for mass in masses]
        # masses = np.logspace(15, 16.5, 100)
        # images = [newPlotter.apply_lensing(BH(mass, 10, 522, 424)) for mass in masses]
        animation(masses, images, 100, output=file.replace('.fits', '_lensing1_1'))

        # # 2. change distance
        # dists = np.linspace(1, 15, 100)
        # images = [newPlotter.apply_lensing(BH(1e15, dist, 300, 300)) for dist in dists]
        # animation(dists, images, 100, output=file.replace('.fits', '_lensing2'))

        # 3. change position
        # xs = np.linspace(100, 500, 201, dtype=int)
        # images = [newPlotter.apply_lensing(BH(1e15, 2, x, 260)) for x in xs]
        # animation(xs, images, 100, output=file.replace('.fits', '_lensing3'))