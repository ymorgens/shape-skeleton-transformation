# **Shape Skeleton transformations**


The following code shows how to create shape variants by manipulating a shape's skeletal representation.  The skeletal representation decomposes the shape into parts.  Here we vary the relationship of the lengths, orientations, widths, and positions of all these parts to create new shape variants.  

Such new variants may be useful in producing stimuli for a variety of psychophysical tasks.  As an example, the code can be used to create classes of objects with similar/dissimilar statistics (e.g., by choosing the distribution of parameters of how these parts are sampled).

## **Requirements**
MATLAB

## **Example**

As an example (See **demo_transformBaseShapeRandom.m**), the novel object and the subparts of its skeleton representation (_right_) are transformed (_left_) to produce a new variant.  Blurring the shape (_lower left_) can remove discontinuities that arise due to the transformation.

<p align="center">
  <img width="460" src="https://github.com/ymorgens/shape-skeleton-transformation/blob/master/skeleton_transformation/demoexample.png">
</p>

## **Reference**
***
If you use this code, please cite the following:

(1) **Morgenstern, Y., Schmidt, F., & Fleming, R. W. (2019). One-shot categorization of novel object classes in humans. Vision Research, 165, 98-108.**

This demo relies on the shape skeleton implementation by Feldman and Singh (2006).  The ShapeToolbox implementation  was downloaded at http://ruccs.rutgers.edu/images/ShapeToolbox1.0.zip

(2) **Feldman, J., & Singh, M. (2006). Bayesian estimation of the shape skeleton. Proceedings of the National Academy of Sciences, 103(47), 18014-18019.**

For questions, please contact Yaniv Morgenstern([Yaniv.Morgenstern@psychol.uni-giessen.de](mailto:Yaniv.Morgenstern@psychol.uni-giessen.de)).

