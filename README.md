BunnyKiller Transient Renderer
==============================

BunnyKiller Transient Renderer is a small physically based renderer implementing  
the rendering algorithms developed for transient rendering described in:  

> "A Framework for Transient Rendering"  
> Adrian Jarabo, Julio Marco, Adolfo Muñoz, Raul Buisan, Wojciech Jarosz, Diego Gutierrez  
> ACM Transactions on Graphics, Vol.33(6), SIGGRAPH Asia 2014  

-----

Based on the paper above, we add some new features for transient rendering of moving objects. 

Overview
--------

The code compiles to a command line program that can render HDR images of scenes  
stored as obj files. The format of the command line arguments is explained in  
INPUT.md.

* The code uses the CImg library for output HDR images. Thus, while no  
  additional linking or dll is needed for output HDR format, other formats can  
  be used by using additional dlls.
* The code is designed with educational purposes, and therefore is not optimized  
  for speed, but to be the most comprehensive and extendible as possible.
* The original [project page is here](http://giga.cps.unizar.es/~ajarabo/pubs/transientSIGA14/code).

## 代码说明

类```RenderEngine```为实现渲染的入口类，会调用以下函数执行渲染：

```c++
RenderEngine->render(const char* name, World, Integrator, Camera, Film, Sampler)
```

下面分别说明参数中的类：

1. ```World```存储整个场景中的模型、光源、材质等信息，实现光线与环境求交的功能；
2. ```Integrator```积分器类，计算渲染方程中的积分。计算过程中对积分进行蒙特卡洛采样，并返回采样后计算得到的所有样本值；
3. ```Camera```相机类，存储相机参数；
4. ```Film```用于将积分器得到的样本，按照一定规则累积到视频中某一帧的某一位置上，最终生成HDR图片；
5. ```Sampler```采样器类，由于整个视频是一个$T\times H\times W \times D$的结构（时间、长、宽、颜色通道），采样器对其中若干个维度采样（比如像素位置x,y），使用积分器在该维度确定情况下采集剩余维度的样本。

