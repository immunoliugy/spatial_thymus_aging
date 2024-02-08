# aging_mouse_thymus

Provided here are the code and functions that we used to create each figure in our manuscript. All data are available on https://singlecell.broadinstitute.org/single_cell/study/SCP2424

## Quickstart
### Install
- You can download Jupyter, here: https://jupyter.org/install. Instructions on running Jupyter notebooks are available here: https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/execute.html.
- Download all of the data here: https://singlecell.broadinstitute.org/single_cell/study/SCP2424/
- Install packages by running the following:
```shell
conda create --name <env> --file requirements.txt
```
- Download the jupyter notebook into the same directory as the data and run each block to re-create all the figures.

## Additional files
- The MiXCR txt file shows an example of the commands used for MiXCR.
- The R code shows how unsupervised clustering and RCTD were performed on the data.
- save_trac_and_trbc2_constant.py shows how the TCR constant UMIs were pulled from the Slide-seq data.
