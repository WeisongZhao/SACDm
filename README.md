
[![website](https://img.shields.io/badge/website-up-green.svg)](https://weisongzhao.github.io/SACDm/)
[![paper](https://img.shields.io/badge/paper-nature%20photon.-black.svg)](https://www.nature.com)
[![Github commit](https://img.shields.io/github/last-commit/WeisongZhao/SACDm)](https://github.com/WeisongZhao/SACDm/)
[![License](https://img.shields.io/github/license/WeisongZhao/SACDm)](https://github.com/WeisongZhao/SACDm/blob/master/LICENSE.txt)<br>
[![Twitter](https://img.shields.io/twitter/follow/weisong_zhao?label=weisong)](https://twitter.com/weisong_zhao/status/1370308101690118146)
[![GitHub watchers](https://img.shields.io/github/watchers/WeisongZhao/SACDm?style=social)](https://github.com/WeisongZhao/SACDm/) 
[![GitHub stars](https://img.shields.io/github/stars/WeisongZhao/SACDm?style=social)](https://github.com/WeisongZhao/SACDm/) 
[![GitHub forks](https://img.shields.io/github/forks/WeisongZhao/SACDm?style=social)](https://github.com/WeisongZhao/SACDm/)

<p>
<h1 align="center">SACD<font color="#b07219">m</font></h1>
<h5 align="center">Fater super-resolution fluctuation imaging: SACD reconstruction with MATLAB.</h5>
<h6 align="right">v0.2.0</h6>
</p>
<br>

<!-- <p>
<img src='./imgs/splash.png' align="left" width=170>
<p> -->

<p>
<img src='./imgs/MATLAB.jpg' align="left" width=120>
</p>

This repository is for SACD reconstruction, and it will be in continued development. It is distributed as accompanying software for publication: [Weisong Zhao et al. Enhancing detectable fluorescence fluctuation for high-throughput and four-dimensional live-cell super-resolution imaging, XXX (2022)](https://www.nature.com). Please cite SACD in your publications, if it helps your research.
<br>
<br>
<br>

<!-- <p>
<img src='./imgs/imagej-128.png' align="right" width=50>
</p>
<br>

[Portal]() to the plugin. -->

## SACD reconstruction

<p align='center'>
<img src='./imgs/Fig1.png' align="center" width=900>
</p>

## Instruction

- The SACD reconstruction requires resolution-related parameter to execute deconvolution, you can give it with `objective-NA`; `wavelength (nm)`; and `pixel-size (nm)`, or just provide `resolution` and `pixel-size`, or feed it with your `own-PSF`. Here are 3 examples:
```python
SRimg = SACDm(imgstack, 'pixel', 65, 'NA', 1.3, 'wavelength', 561);
SRimg = SACDm(imgstack, 'pixel', 65, 'resolution', 250);
SRimg = SACDm(imgstack, 'psf', ownpsf);
```

- Please try help to get the API.
```python
addpath(genpath('SACDm')); 
help SACDm
```

- Regarding the SACD SR frame visualization, it can be scaled with a gamma correction according to the bSOFI setting.
```python
background = 0.02; order = 2;
SRimg2vis = SRimg.^0.5;
SRimg2vis(SRimg2vis < order * background * max(SRimg2vis(:))) = 0;
figure(2);imshow(SRimg2vis, [], 'colormap', hot)
```

Two demos can also be found at the [SACDj release v1.1.3](https://github.com/WeisongZhao/SACDj/releases/tag/v1.1.3).


## Declaration
This repository contains the MATLAB source code for <b>SACD</b> .  

If you are not a MATLAB user, you can have a try on the imagej version of SACD: [SACDj](https://github.com/WeisongZhao/SACDj).

<p>
<img src='./imgs/imagej-128.png' align="right" width=50>
</p>
<br>
<br>

## Version
- v0.2.0 Sparse-SACD reconstruction core
- v0.1.0 SACD reconstruction core

## Related links: 
- ImageJ plug-in version of SACD: [SACDj](https://github.com/WeisongZhao/SACDj)
- **Some fancy results and comparisons:** [my website](https://weisongzhao.github.io/MyWeb2/portfolio-4-col.html)
- **Preprint:** [Weisong Zhao et al. Enhancing detectable fluorescence fluctuation for high-throughput and four-dimensional live-cell super-resolution imaging, bioRxiv (2022).](https://doi.org/10.1101/2022.12.12.520072)
- **Reference:** [COMING SOON](#)

<details>
<summary><b>Plans</b></summary>

- Full FRC assisted SACD;
- Full 3D-SACD;
- GPU acceleration.
</details>

## Open source [SACDm](https://github.com/WeisongZhao/SACDm)
This software and corresponding methods can only be used for **non-commercial** use, and they are under Open Data Commons Open Database License v1.0.