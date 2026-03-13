import matplotlib.pyplot as plt
from pathlib import Path

from BLPlot.plotter import (
    Plotter,
    iter_datasets_with_runs,
    load_dataset_metric,
    make_box_figure,
)


class PlotEPR(Plotter):
    """
    Plotter that produces one early precision box plot per dataset.

    For each dataset, reads EarlyPrecision.csv and draws a box plot with
    one box per algorithm showing the distribution across runs. A dashed
    reference line marks the random-predictor EPR baseline. Each dataset is
    written as both a PDF and PNG under an EPR/ subdirectory of the output
    directory, named <dataset_label>-EPR.pdf and <dataset_label>-EPR.png.
    """

    def __call__(self, config: dict, output_dir: Path, root: Path) -> None:
        """
        Generate per-dataset EPR box plots, writing each to its own PDF and PNG.

        Creates output_dir/EPR/ and writes <dataset_label>-EPR.pdf and
        <dataset_label>-EPR.png for each enabled dataset. Datasets that
        produce no figure (missing data) are skipped.

        Parameters
        ----------
        config : dict
            Parsed YAML configuration.
        output_dir : Path
            Parent directory; an EPR/ subdirectory is created inside it.
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

        epr_dir = output_dir / 'EPR'
        epr_dir.mkdir(parents=True, exist_ok=True)

        for _, dataset_label, dataset_path, gt_path, _ in iter_datasets_with_runs(config, root):
            fig = make_box_figure(
                load_dataset_metric(dataset_path, 'EarlyPrecision.csv'),
                f'Early Precision Ratio — {dataset_label}', 'EPR',
                rand_value=1.0,
                ylim=None,
            )
            if fig is None:
                continue
            safe_label = dataset_label.replace('/', '-')
            stem = epr_dir / f'{safe_label}-EPR'
            fig.savefig(stem.with_suffix('.pdf'))
            fig.savefig(stem.with_suffix('.png'))
            plt.close(fig)
            print(f"Saved {stem}.pdf and .png")
