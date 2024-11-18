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

```python
class BH:
    def __init__(self, mass: float, distance: float, x: float, y: float):
        self.mass = mass            # Solar masses
        self.distance = distance    # Mpc
        self.x = x                  # x-coordinate(pixel)
        self.y = y                  # y-coordinate(pixel)
    
    def __str__(self) -> str:
        return f'BH(mass={self.mass}, distance={self.distance}, x={self.x}, y={self.y})'
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
