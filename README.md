# Partial Multi view clutering
Matlab Codes for Partial Multi view clutering methods such as MIC, IMG, GPMVC, MKKIK, PVC and my PMVC-ASC(PMLAN).

Repository for the paper "Partial Multi-view Clustering via Auto-Weighting Similarity Completion"(https://link.springer.com/chapter/10.1007/978-3-319-97909-0_23).

## Abstract

With the development of data collection techniques, multi-view clustering (MVC) becomes an emerging research direction to improve the clustering performance. However, most MVC methods assume that the objects are observed on all the views. As a result, existing MVC methods may not achieve satisfactory performance when some views are incomplete. In this paper, we propose a new MVC method, called as partial multi-view clustering via auto-weighting similarity completion (PMVC-ASC). The major contribution lies in jointly learning the consensus similarity matrix, exploring the complementary information among multiple distinct feature sets, quantifying the contribution of each view and splitting the similarity graph into several informative submatrices, each submatrix corresponding to one cluster. The learning process can be modeled via a joint minimization problem, and the corresponding optimization algorithm is given. A series of experiments are conducted on real-world datasets to demonstrate the superiority of PMVC-ASC by comparing with the state-of-the-art methods.

## Citation

If you find this project useful in your research, please consider cite:

@inproceedings{min2018partial,
  title={Partial Multi-view Clustering via Auto-Weighting Similarity Completion},
  author={Min, Chen and Cheng, Miaomiao and Yu, Jian and Jing, Liping},
  booktitle={Chinese Conference on Biometric Recognition},
  pages={214--222},
  year={2018},
  organization={Springer}
}

## License

Our code is released under the Apache 2.0 license.
