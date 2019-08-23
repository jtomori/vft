VFX Fractal Toolkit
==========================
![VFX Fractal Toolkit banner image](img/vft_cover.jpg)
*Set of tools for generating fractal and generative art.*

<br>

## About
This is my graduation project: **VFX Fractal Toolkit** (VFT), which I developed at *Filmakademie Baden-Württemberg* while studying [Technical Directing](https://animationsinstitut.de/de/studium/animation/technical-director/informationen/).

It contains tools written in *OpenCL, OSL, Blink, Python, VEX and JavaScript* intended to be used in *Houdini, Arnold, Nuke or a web browser*.

The code is in prototyping stage and many features are experimental. It is **not production ready** and most parts of it need refactoring.

Here are some animations produced with it:

* [![Volumetric fractals](img/volumes.jpg)](https://www.youtube.com/watch?v=E8n6chN2Txw)
* [![Dynamical systems](img/particles.jpg)](https://www.youtube.com/watch?v=_gdApm_QPjs)
* [![2D fractals](img/2d.jpg)](https://www.youtube.com/watch?v=__8gaEv5GAs)

I had a chance to present progress of VFT at two **FMX** conferences (2018, 2019), you can find the recordings here:

* [![FMX 2019 recording](img/fmx_19.jpg)](https://youtu.be/n-m00N7TYYM?t=2452)
* [![FMX 2018 recording](img/fmx_18.jpg)](https://youtu.be/SNa18n5d8UY?t=1m26s)

It was also featured in **Posters Preview: SIGGRAPH 2019** video:
* [![Siggraph 2019 posters preview](img/sig_19.jpg)](https://youtu.be/aRmfaEBLNmw)

This project was presented at **The 15th ACM SIGGRAPH European Conference on Visual Media Production** [(CVMP 2018)](https://www.cvmp-conference.org/2018/programme/) conference: [fast-forward](https://www.youtube.com/watch?v=_CI8GFDmKZQ), [paper](https://animationsinstitut.de/fileadmin/user_upload/files_forschung/pdf/Publications/18_cvmp_vft_juraj_tomori_paper.pdf), [poster](https://animationsinstitut.de/fileadmin/user_upload/files_forschung/pdf/Publications/18_cvmp_vft_juraj_tomori_poster.png).

It was also presented in [posters session](https://s2019.siggraph.org/presentation/?sess=sess175&id=pos_114#038;id=pos_114) at **SIGGRAPH 2019** in Los Angeles. You can find the abstract [here](https://dl.acm.org/citation.cfm?id=3306214.3338543).

You can also cite my work:
```
@inproceedings{Tomori:2019:VFT:3306214.3338543,
 author = {Tomori, Juraj},
 title = {VFX Fractal Toolkit: Integrating Fractals into VFX Pipeline},
 booktitle = {ACM SIGGRAPH 2019 Posters},
 series = {SIGGRAPH '19},
 year = {2019},
 isbn = {978-1-4503-6314-3},
 location = {Los Angeles, California},
 pages = {97:1--97:2},
 articleno = {97},
 numpages = {2},
 url = {http://doi.acm.org/10.1145/3306214.3338543},
 doi = {10.1145/3306214.3338543},
 acmid = {3338543},
 publisher = {ACM},
 address = {New York, NY, USA},
 keywords = {distance fields, fractals, pipeline, ray marching, vfx},
} 
```

You can find **comparison of various techniques** (visual quality vs performance) [here](comparison.md).

<br>

## Resources
* [Mandelbulber2 source code repository](https://github.com/buddhi1980/mandelbulber2/)
* [Mandelbulb3D source code repository](https://github.com/thargor6/mb3d)
* [Capturing the infinite universe in "Lucy": fractal rendering in film production](https://dl.acm.org/citation.cfm?id=2614166)
* [The fractal nature of Guardians of the Galaxy Vol. 2](https://www.fxguide.com/featured/the-fractal-nature-of-guardians-of-the-galaxy-vol-2/)

<br>

## Thanks
* [Íñigo Quílez](http://www.iquilezles.org/www/index.htm) - great articles on raymarching, fractals, orbit traps, SDFs...
* Krzysztof Marczak - lead Mandelbulber2 developer, supporting via emails
* [Mikael Hvidtfeldt Christensen](http://blog.hvidtfeldts.net/) - great articles on raymarching, fractals, generative art
* [Dom Penfold](http://woo4.me/) - blog with useful articles