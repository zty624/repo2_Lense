# 天体物理第二次程序作业

题目要求：实现一个简单的引力透镜的程序模拟效果。

## 数据获取和处理

从 NED 中下载的数据是 `FITS` 格式的，我们需要使用 `astropy` 库来处理这些数据。

```python
from astropy.io improt fits
import os
DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
(_, _, DATA_FILES) = next(os.walk(DATA_DIR))
def read(file: str) -> tuple[2]:
    with fits.open(file) as h:
        header, data = h[0].header, h[0].data
    return (header, data)
```

使用 `plt.imshow` 绘制图像。
```python
class Plotter:
    def __init__(self, data: tuple):
        self.header, self.data = data

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
```

得到原图像的灰度图：
![original](https://github.com/zty624/repo2_Lense/blob/main/src/data/NGC_1300_I_IIIaJ_dss1.png)

## 引力透镜的代码实现

### 单位换算

考虑到照片是以像素为单位，而实际上我们只能得到像素和观测角度的关系，所以我们需要将像素转换为角度，然后在平面近似下进行光线追踪。
在 `FITS` 的头文件中搜索关于像素和角度的分辨率的数据，换算得到每个 `pixel` 和角度的比例尺。

```python
class Plotter:
    def __init__(self, data: tuple):
        self.header, self.data = data
        xpixelsz = self.metadata.loc[self.metadata['Key'] == 'XPIXELSZ'].iloc[0]['Value']
        ypixelsz = self.metadata.loc[self.metadata['Key'] == 'YPIXELSZ'].iloc[0]['Value']
        pltscale = self.metadata.loc[self.metadata['Key'] == 'PLTSCALE'].iloc[0]['Value']
        self.xscale = xpixelsz * pltscale / 1000 * np.pi / 180 / 3600   # rad/px
        self.yscale = ypixelsz * pltscale / 1000 * np.pi / 180 / 3600   # rad/px

    @property
    def metadata(self):
        return pd.DataFrame(self.header.items(), columns=['Key', 'Value'])
```

### 动画实现

使用 `matplotlib.animation` 包绘制 gif 动图。

```python
import matplotlib.animation as ani
def animation(images: list, output: str = None):
    fig, ax = plt.subplots(figsize=(12, 12))
    ax.axis('off')
    im = ax.imshow(images[0], cmap='gray')
    def update(i):
        im.set_data(images[i])
        ax.set_title(f'Mass: {mass_list[i]:.2e} kg')
        return im,
    anime = ani.FuncAnimation(fig, update, frames=len(images), interval=200)
    if output:
        anime.save(os.path.join(DATA_DIR, output + ".gif"), writer='imagemagick')
```

### 引力透镜变换

本项目作如下近似：

1. 原图像在天球上所张的立体角远小于 1，以至于在此范围内可以使用平面近似。
2. 黑洞距离地球足够远，使得旁轴条件成立，可以使用简单的方程。

本题建模如下：

1. 采用简单的反向光线追迹法，假想望远镜发出了射向原目标范围 $(-\theta_0 \sim \theta_0, -\theta_0 \sim \theta_0)$ 的光线集。
2. 对光线集施加一个黑洞带来的额外的偏折，使得原来朝向 $(\theta_x, \theta_y)$ 的光线最终射到了 $(\theta_x', \theta_y')$ 的背景上。
3. 将 $(\theta_x', \theta_y')$ 上的亮度投影到 $(\theta_x, \theta_y)$ 的坐标上。

```python
class Plotter:
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
```

## 输出结果

此处缺失 gif