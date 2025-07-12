# MW-PFCM-Brain-MRI
An Adaptive Multi-Weighted Possibilistic Fuzzy C-Means Clustering Approach for Brain Tumor Segmentation in Magnetic Resonance Imaging

The possibilistic fuzzy c-means (PFCM) method is a well-known approach to practical image segmentation. However, it still has some limitations, including 1) considering equal importance to all image features, 2) sensitivity to initialization, and 3) ignoring the neighborhood information of image pixels in the objective function. To alleviate the mentioned challenges, this paper proposes a novel adaptive multi-weighted PFCM image segmentation approach. The proposed possibilistic fuzzy objective function is composed of 1) a feature weighting factor to distinguish between the importance of different features, 2) a cluster weighting schema to mitigate the initialization sensitivity, 3) a spatial constraint term containing local information of pixels to reduce the algorithm's sensitivity to noise. In addition, a non-Euclidean distance metric based on a kernel is utilized to improve robustness to feature scale differences and better capture nonlinear structures within the data. Furthermore, the curvelet transform is employed in this research to enhance the edges in brain Magnetic Resonance Imaging (MRI) images, their contrast, and brightness intensity. Moreover, an efficient combination of image features increases the segmentation accuracy. Finally, mathematical analyses are provided to obtain the required updating functions. In summary, the artificial intelligence-powered segmentation method integrates adaptive weighting and spatial constraints to better process and interpret complex brain images and, as a result, develop efficient computer-aided diagnosis (CAD) systems. The proposed method's performance is evaluated and compared with state-of-the-art methods. The results show that the proposed method outperforms the competitors.


# Overview of MW-PFCM-Brain-MRI:
<img width="630" height="215" alt="Untitled" src="https://github.com/user-attachments/assets/23ed51b3-a276-43e1-8ee9-48011f6ea65b" />


# Comment:
The repository file includes the MATLAB implementation of the mw_pfcm_brain_mri algorithm.

Comments are written for all steps of the algorithm for better understanding the code. Also, a demo is implemented for ease of running, which runs by importing the data and other necessary algorithm parameters.

# Used Dataset:
To evaluate the proposed method, the benchmark dataset has been used. There are some images of this dataset in the uploaded file (Dataset folder). Whole dataset is available on: https://www.kaggle.com/datasets/mateuszbuda/lgg-mri-segmentation.

## Condition and terms to use any sources of this project (Codes, Datasets, etc.):

1) Please cite the following papers:

[1] M. Hashemzadeh, A. Golzari Oskouei, and N. Farajzadeh, "New fuzzy C-means clustering method based on feature-weight and cluster-weight learning," Applied Soft Computing, vol. 78, pp. 324-345, 2019/05/01/ 2019, doi: https://doi.org/10.1016/j.asoc.2019.02.038.

[2] A. Golzari Oskouei, M. Hashemzadeh, B. Asheghi, and M. A. Balafar, "CGFFCM: Cluster-weight and Group-local Feature-weight learning in Fuzzy C-Means clustering algorithm for color image segmentation," Applied Soft Computing, vol. 113, p. 108005, 2021/12/01/ 2021, doi: https://doi.org/10.1016/j.asoc.2021.108005.

[3] A. Golzari Oskouei and M. Hashemzadeh, "CGFFCM: A color image segmentation method based on cluster-weight and feature-weight learning," Software Impacts, vol. 11, p. 100228, 2022/02/01/ 2022, doi: https://doi.org/10.1016/j.simpa.2022.100228.

2) Please do not distribute the database or source codes to others without the authorization from Dr. Mahdi Hashemzadeh.

Authors’ Emails: arezounajafi[at]azaruniv.ac.ir (A.N. Moghadam), nasseraghazadeh[at]iyte.edu.tr (N. Aghazadeh), hashemzadeh[at]azaruniv.ac.ir (M. Hashemzadeh), and a.golzari[at]tabrizu.ac.ir (A. G. Oskouei)
