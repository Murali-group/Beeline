from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from BLPlot.plotter import (
    Plotter,
    iter_datasets_with_runs,
    load_dataset_metric,
    make_box_figure,
    random_classifier_baseline,
)


class PlotEPR(Plotter):
    """
    Plotter that produces one early precision box plot per dataset.

    For each dataset, reads EarlyPrecision.csv and draws a box plot with
    one box per algorithm showing the distribution across runs. A dashed
    reference line marks the random-predictor EPR baseline. All pages are
    written to a single EPR.pdf in the output directory.
    """

    def __call__(self, config: dict, output_dir: Path, root: Path) -> None:
        """
        Generate per-dataset EPR box plots and write all pages to EPR.pdf.

        Parameters
        ----------
        config : dict
            Parsed YAML configuration.
        output_dir : Path
            Directory where EPR.pdf is written.
        root : Path
            Working directory from which config paths are resolved.

        Returns
        -------
        None
        """
        if not isinstance(config, dict):
            raise TypeError(f"config must be dict, got {type(config)}")
        if not isinstance(output_dir, Path):
            raise TypeError(f"output_dir must be Path, got {type(output_dir)}")
        if not isinstance(root, Path):
            raise TypeError(f"root must be Path, got {type(root)}")

        out_path = output_dir / 'EPR.pdf'
        pages_written = 0

        with PdfPages(out_path) as pdf:
            for dataset_id, dataset_path, gt_path, _ in iter_datasets_with_runs(config, root):
                baseline = (
                    random_classifier_baseline(gt_path)
                    if gt_path.exists() else None
                )
                values = load_dataset_metric(dataset_path, 'EarlyPrecision.csv')
                fig = make_box_figure(
                    values, f'Early Precision — {dataset_id}', 'Early Precision',
                    rand_value=baseline,
                )

                if fig is None:
                    continue
                pdf.savefig(fig)
                plt.close(fig)
                pages_written += 1

        if pages_written:
            print(f"Saved {pages_written} plot(s) to {out_path}")
        else:
            print(f"No EPR data found; {out_path} not written.")
