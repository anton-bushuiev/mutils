import numpy as np
import pandas as pd

from mutils.matrices import ScoringMatrix, ScoringMatrixCollection


def test_collection_average():
    matrices = ScoringMatrixCollection([
        ScoringMatrix(
            pd.DataFrame(
                [
                    [1, 2, 3],
                    [4, 5, 6],
                ],
                columns=['A', 'B', 'C']
            )
        ),
        ScoringMatrix(
            pd.DataFrame(
                [
                    [7, 8, 9],
                    [np.nan, 11, 12],
                ],
                columns=['A', 'B', 'C']
            )
        )
    ])

    expected_unit_weights = ScoringMatrix(
        pd.DataFrame(
            [
                [4., 5., 6.],
                [4., 8., 9.],
            ],
            columns=['A', 'B', 'C']
        )
    )
    assert (expected_unit_weights == matrices.average())

    expected_weighted = ScoringMatrix(
        pd.DataFrame(
            [
                [5., 6., 7.],
                [4., 9., 10.],
            ],
            columns=['A', 'B', 'C']
        )
    )
    matrices.matrices[1].weight = 2.
    assert (expected_weighted == matrices.average())

    matrices.matrices[1].weight = 1.
    assert (expected_unit_weights == matrices.average())
