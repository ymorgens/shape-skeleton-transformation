# **Shape Skeleton transformations**


The following code shows how to create shape variants by manipulating a shape's skeletal representation.  The skeletal representation decomposes the shape into parts.  Here we vary the relationship of the lengths, orientations, widths, and positions of all these parts to create new shape variants.  

Such new variants may be useful in producing stimuli for a variety of psychophysical tasks.  As an example, the code can be used to create classes of objects with similar/dissimilar statistics (e.g., by choosing the distribution of parameters of how these parts are sampled).

<p align="center">
  <img width="460" src="https://github.com/ymorgens/shape-skeleton-transformation/blob/master/skeleton_transformation/readme image.png">
</p>

## **Requirements**
MATLAB

## **Example**

As an example (See **demo_transformBaseShapeRandom.m**), the novel object and the subparts of its skeleton representation (_right_) are transformed (_left_) to produce a new variant.  Blurring the shape (_lower left_) can remove discontinuities that arise due to the transformation.

<p align="center">
  <img width="460" src="https://github.com/ymorgens/shape-skeleton-transformation/blob/master/skeleton_transformation/demoexample.png">
</p>

The code may be useful in producing novel object classes for psychophyisical experiments.  For example, one can create different classes of objects by manipulating the distribution that the skeletal parameters are drawn from (see **demo_ShapeEnvs.m**)
<p align="center">
  <img width="460" src="https://github.com/ymorgens/shape-skeleton-transformation/blob/master/skeleton_transformation/demo_ShapeEnvs1.png">
</p>

## **References**
***
If you use this code, please cite the following:

(1) **Morgenstern, Y., Schmidt, F., & Fleming, R. W. (2019). [One-shot categorization of novel object classes in humans](https://www.sciencedirect.com/science/article/abs/pii/S0042698919301749). Vision Research, 165, 98-108.**

(2) **Morgenstern, Y., Schmidt, F., & Fleming, R. W. (2020). [A dataset for evaluating one-shot categorization of novel object classes](https://www.sciencedirect.com/science/article/pii/S2352340920301967). Data in brief, 29, 105302.** [[Data]](https://zenodo.org/record/3433278#.X1ZEgdMzaJQ)

This demo relies on the shape skeleton implementation by Feldman and Singh (2006).  The ShapeToolbox implementation  was downloaded at http://ruccs.rutgers.edu/images/ShapeToolbox1.0.zip

(3) **Feldman, J., & Singh, M. (2006). [Bayesian estimation of the shape skeleton. Proceedings of the National Academy of Sciences](https://www.pnas.org/content/pnas/103/47/18014.full.pdf), 103(47), 18014-18019.**

## **Related work**
***
For related work, please visit the sites of [Roland Fleming](http://fleming.oerloeg.com/),[Filipp Schmidt](http://www.allpsych.uni-giessen.de/filipp/), and [Yaniv Morgenstern](https://sites.google.com/view/yanivmorgenstern) 


## **Questions**
***
For questions, please contact Yaniv Morgenstern ([Yaniv.Morgenstern@psychol.uni-giessen.de](mailto:Yaniv.Morgenstern@psychol.uni-giessen.de)).

