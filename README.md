# **spatial_thymus_aging**
Provided here are the code and functions that we used to create each figure in our manuscript. All data are available on [Single Cell Portal SCP2424](https://singlecell.broadinstitute.org/single_cell/study/SCP2424).

---

## **Quickstart with Code**
### **Install**
1. **Install Jupyter Notebook**  
   - You can download Jupyter [here](https://jupyter.org/install).  
   - Instructions on running Jupyter notebooks are available [here](https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/execute.html).  

2. **Download the Data**  
   - The dataset can be found at: [Single Cell Portal SCP2424](https://singlecell.broadinstitute.org/single_cell/study/SCP2424/)  
   - Download all required files into a directory.

3. **Install Required Packages**  
   - Create a Conda environment and install dependencies:
   ```shell
   conda create --name <env> --file requirements.txt
   conda activate <env>
   ```

---

## **Quickstart with Code Ocean**
This Code Capsule is hosted on **Code Ocean**, a platform recommended for **reproducible code**.

---

### **Option 1: Run All Code and Output Results (Fast, Recommended)**
This option **runs all scripts** and **generates all figures**, saving them in the `results` folder. It is the best option for verifying figure outputs efficiently.

#### **Steps:**
1. **Set the script to run:**
   - Locate the file **`run`** under the `code` directory.
   - Right-click it and select **"Set as File to Run"**.
   - _(Optional: To run a specific figure, modify this file to comment out other figure commands.)_

2. **Start the Reproducible Run:**
   - Click **"Reproducible Run"**.
   - Main figures (**Figs. 1–5**) will take **~50 minutes** to complete.
   - You may exit the window and return later to view the results.

3. **Review the results:**
   - Open the `results` folder to access:
     - All **generated figures**.
     - The executed **Jupyter Notebook** saved as an **HTML file** (with outputs for easy review).
   - **Download results** if needed before running supplementary figures.

4. **Run Supplementary Figures (Optional):**
   - Repeat **steps 1–3**, but this time, right click and set **`run_supp`** as the file to run.
   - Supplementary figures will take **~40 minutes** to complete.

---

### **Option 2: Run Each Figure Manually (Slow, for Detailed Code Review)**
This method allows for **step-by-step execution** but may require restarting due to **memory constraints**.

#### **Steps:**
1. **Launch JupyterLab:**
   - Under **"Reproducible Run"**, click **"Launch a Cloud Workstation"** and select the **JupyterLab** icon.

2. **Wait for Capsule Initialization:**
   - The environment includes all necessary data, code, and dependencies and may take a few minutes to load.

3. **Set Up the Environment:**
   - Once loaded, execute the code in the **“Set up environment”** section.
   - _(This step is required before running any analysis.)_

4. **Run Figure-Specific Code:**
   - Open the notebook for the desired figure.
   - Go to **Kernel > Restart Kernel and Run All Cells…**
   - _(Note: If the kernel dies due to memory limitations, re-run the "Set up environment" code before executing any figure-specific code again.)_
   - You do not need to run previous figure scripts to execute a specific panel.
     _(For example, you can run **Fig. 4I** independently without running **Figs. 4F-H**.)_
