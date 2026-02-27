import sys
import multiprocessing

# Force the mapping of importlib_metadata to importlib.metadata
try:
    import importlib.metadata
except ImportError:
    import importlib_metadata
    sys.modules['importlib.metadata'] = importlib_metadata

# Apply patch: Skip deleting temporary folders
def apply_patch():
    import os
    import shutil

    original_remove_dir = None

    def patched_remove_dir(data_path, adata):
        """
        The patched version of the "remove_dir" function, skipping the deletion step
        Only creates folders, without deleting existing folders
        """
        dir = os.path.split(data_path)[0]
        filename = os.path.splitext(os.path.basename(data_path))[0]
        NEW_DIR = os.path.join(dir, filename)
        adata.uns['temp'] = NEW_DIR

        # Only creates folders, without deleting existing folders
        if not os.path.exists(NEW_DIR):
            os.mkdir(NEW_DIR)

        log_file = os.path.join(NEW_DIR, "PATCH_INFO.txt")
        with open(log_file, 'w') as f:
            f.write(f"The temporary folder has been created but not deleted. Please manually delete this folder：{NEW_DIR}\n")
            f.write(f"Creation time：{os.path.getctime(NEW_DIR) if os.path.exists(NEW_DIR) else 'N/A'}\n")

        return NEW_DIR

    try:
        import unitvelo.utils as utils_module
        # Save the original function
        original_remove_dir = utils_module.remove_dir
        # Apply the patch
        utils_module.remove_dir = patched_remove_dir
        print("patch applied: skipping the deletion step")
    except Exception as e:
        print(f"applying patch failed：{e}")

if __name__ == '__main__':
    multiprocessing.freeze_support()
    apply_patch()

import unitvelo as utv
import scvelo as scv
import numpy as np
import matplotlib.pyplot as plt

def main():
    adata = scv.read("../../data/scPDAC_with_su_normalized.h5ad")

    velo = utv.config.Configuration()
    velo.R2_ADJUST = True
    velo.FIT_OPTION = '1'  # Unified-time
    velo.GPU = 0

    adata = utv.run_model(adata, 'condition', config_file=velo)
    adata.write_h5ad("../../results/UniTVelo/adata/adata_with_UniTVelo.h5ad")

if __name__ == '__main__':
    main()