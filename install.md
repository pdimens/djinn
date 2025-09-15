---
label: Installing
description: Installing Djinn
icon: desktop-download
order: 99
---

## conda
!!!warning
**NOT YET AVAILABLE**
!!!

Djinn can be installed with conda/mamba:
```bash
conda install -c bioconda -c conda-forge bioconda::djinn
```


## pip
Alternatively, if conda isn't available on your system, you can
install the python package manually. You will also need to make sure `samtools` is installed on your system.
### 1. get the latest release
This example uses version `1.0`, but the version you want may differ
```bash
DJINN_VERSION="1.0"
https://github.com/pdimens/djinn/releases/download/${DJINN_VERSION}/djinn.${DJINN_VERSION}.tar.gz
```

### 2. decompress it
```bash
tar xvfz djinn.${DJINN_VERSION}.tar.gz
cd djinn.${DJINN_VERSION}
```

### 3. install with `pip`
```bash
python -m pip install .
```

```bash in the event the installation wont work, try this
python3 -m pip install --upgrade build && python3 -m build && \
    pip install dist/*.whl
```
