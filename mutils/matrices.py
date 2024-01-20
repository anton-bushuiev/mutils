"""
TODO Fix the synchronization of list and dict
```
for matrix in matrices:
    if matrix.name == 'MutCompute2 (1BUI)':
        display(matrix.df)
display(matrices['MutCompute2 (1BUI)'].df)
```
"""

from __future__ import annotations

import copy
import operator
from pathlib import Path
from collections import defaultdict
from typing import Iterable, Collection, Callable, List, Any, Union, Dict, \
    Literal

import numpy as np
import pandas as pd
import torch
import matplotlib.pyplot as plt
import seaborn as sns
import plotly
from scipy import stats
from matplotlib.patches import Rectangle

from mutils.proteins import AMINO_ACID_CODES_1, MISSING_AMINO_ACID_CODE


class Matrix:
    def __init__(self):
        raise NotImplementedError()


class ScoringMatrix:
    SCORE_TITLE = 'Score'

    def __init__(
        self,
        df: pd.DataFrame,
        wts: str = None,
        name: str = '',
        weight: float = 1.0,
    ):
        self.df = df
        self.name = name
        self.weight = weight

        self.n_residues = len(self.df)
        self.n_letters = len(self.df.columns)

        if wts is None:
            wts = self.n_residues * MISSING_AMINO_ACID_CODE
        self.wts = wts

    @classmethod
    def from_excel(
        cls,
        path: Path,
        wts: str,
        name: str,
        *args,
        **kwargs
    ):
        # Read values from Excel
        df = pd.read_excel(path, sheet_name=name, index_col=0)

        # Init matrix
        matrix = cls(
            df=df,
            wts=wts,
            name=name,
            *args,
            **kwargs
        )
        return matrix

    @classmethod
    def from_template(
        cls,
        df: Collection,
        template_matrix: ScoringMatrix
    ):
        matrix = copy.deepcopy(template_matrix)
        matrix.df = pd.DataFrame(
            df,
            columns=matrix.df.columns,
            index=matrix.df.index
        )
        return matrix

    def to_excel(self):
        raise NotImplementedError()

    def standardize(
        self,
        mean: float = None,
        std: float = None,
        shift: bool = True
    ):
        # Estimate mean if to shift the values
        if not shift:
            mean = 0.0
        elif mean is None:
            mean = self.df.stack().mean()

        # Estimate standard deviation
        if std is None:
            std = self.df.stack().std()

        # Rescale
        self.df = (self.df - mean) / std

    def subtract_wts(self):
        # Subtract wild types
        def substract_wt(row):
            wt = self.wts[row.name]
            if wt == MISSING_AMINO_ACID_CODE:
                return row
            else:
                return row - row[wt]

        self.df = self.df.apply(substract_wt, axis=1)

    def filter_disruptive_mutations(self):
        self.filter_with_mask(self.df >= 0)
        return self

    def filter_identity_mutations(self):
        mask = np.full_like(self.df, fill_value=True).astype(bool)
        mask = pd.DataFrame(mask, columns=self.df.columns)
        for i, wt in enumerate(self.wts):
            if wt != MISSING_AMINO_ACID_CODE:
                mask.loc[i, wt] = False
        self.filter_with_mask(mask)
        return self
    
    def filter_positions(self, positions: Iterable[int], negate: bool = False):
        mask = copy.deepcopy(self.df)
        mask.loc[:, :] = False
        if negate:
            mask.loc[positions, :] = True
        else:
            mask.loc[~mask.index.isin(positions), :] = True

        self.filter_with_mask(mask)
        return self

    def filter_with_mask(self, mask: Union[np.array, pd.DataFrame]):
        mask = pd.DataFrame(mask, columns=self.df.columns)
        self.df = self.df[mask]
        return self
    
    def filtermap(self, func: Callable[Any, bool]):
        self.df = self.df[self.df.applymap(func)]
        return self

    def revert(self):
        self.df = -self.df
        return self

    def replace_zeros(self, value: Any = None):
        # Replace zeros with min observed non-zero value if not specified
        if value is None:
            matrix_values = self.df.stack()
            value = np.min(matrix_values[matrix_values > 0])

        self.df = self.df.replace({0: value})
        return self

    def log(self):
        # Assert input
        if np.sum(self.df.stack() < 0) > 0:
            raise ValueError('Negative values in matrix.')

        # Apply log
        self.df = np.log(self.df)
        return self

    def rank(self):
        raise NotImplementedError()

    def to_normalized_improvements(
        self,
        normalized: bool = False,
        improvements: bool = False,
    ):
        matrix = NormalizedImprovementScoringMatrix(
            df=self.df,
            wts=self.wts,
            name=self.name,
            weight=self.weight,
            normalized=normalized,
            improvements=improvements
        )
        return matrix

    def predict(self, mutations: Iterable[tuple]) -> List:
        # Currently assumes simple additive model (no epistasis)
        predictions = []
        for mp_mutation in mutations:
            sp_predictions = [self.df.loc[pos, mut] for pos, mut in mp_mutation]
            predictions.append(np.mean(sp_predictions))
        return predictions

    def stack(self, score_title=None) -> pd.DataFrame:
        if score_title is None:
            score_title = self.SCORE_TITLE

        return self.df.stack(dropna=False).to_frame(score_title) \
            .reset_index(names=['Position', 'Amino Acid'])

    def plot_interactive(
        self,
        score_title: str = None,
        hover_text: Iterable = None,
        transparent_hover: bool = False
    ):
        # Process arguments
        if score_title is None:
            score_title = self.SCORE_TITLE

        # Stack matrix to triples
        df_stack = self.stack(score_title)
        x = df_stack['Position']
        y = df_stack['Amino Acid']
        z = df_stack[score_title]
        if hover_text is None:
            hover_text = z.astype(str)

        # Init hover titles to mutations
        wts = self.wts
        hover_title = df_stack.apply(
            lambda row: '<b>' + wts[row['Position']]
                        + str(row['Position']) + row['Amino Acid']
                        + '</b>',
            axis=1
        )
        hover_text = hover_title + '<br>' + hover_text

        # Plot
        data = [plotly.graph_objs.Heatmap(
            x=x.astype(str),  # String is crucial for tickson='boundaries'
            y=y,
            z=z,
            colorscale=[
                [0., '#df0101'],
                [0.5, '#f5f6ce'],
                [1., '#31b404']
            ],
            zmid=0,
            hoverinfo='text',
            text=hover_text,
            colorbar=dict(
                title=score_title,
                title_font_size=10,
                thickness=3,
                titleside='right',
                tickfont=dict(
                    size=9
                )
            ),
        )]

        # width = 1100
        # height = int(width / x.nunique()) * y.nunique()
        layout = plotly.graph_objs.Layout(
            autosize=False,
            width=1100,
            height=250,
            margin=dict(l=20, r=20, t=20, b=20),
            xaxis=dict(
                tickangle=-90,
                tickmode='linear',
                tick0=0,
                dtick=1,
                tickfont=dict(
                    size=9
                )
            ),
            xaxis_nticks=x.nunique(),
            yaxis_nticks=y.nunique(),
            yaxis=dict(
                visible=True,
                autorange='reversed',
                tickfont=dict(
                    size=8
                )
            ),
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            title=self.name,
            title_font_size=13,
            title_x=0.5,
        )
        fig = plotly.graph_objs.Figure(data=data, layout=layout)
        fig.update_coloraxes(cmin=z.min(), cmax=z.max())

        # Setup grid
        fig.update_xaxes(
            showgrid=True,
            gridwidth=0.01,  # Does not work
            gridcolor='rgb(230,230,230)',
            tickson='boundaries'
        )
        fig.update_yaxes(
            showgrid=True,
            gridwidth=0.01,  # Does not work
            gridcolor='rgb(240,240,240)',
            tickson='boundaries'
        )

        # Trick to make hovered cell visible
        if transparent_hover:
            fig.update_layout(
                hovermode='x unified',
                hoverlabel=dict(bgcolor='rgba(255,255,255,0.75)',
                font=dict(color='black')),
                legend=dict(title=None)
            )
        fig.update_xaxes(showspikes=True)
        fig.update_yaxes(showspikes=True)

        return fig

    def plot(
        self,
        cmap: Any = sns.cm.rocket,
        wts: bool = False,
        value_label: str = '',
        labels: Dict[tuple, str] = None,
        gaps: bool = False
    ):
        # Init figure for a heatmap with squared cells
        factor = 4
        plt.figure(
            figsize=(self.n_residues // factor, self.n_letters // factor)
        )

        # Plot
        ax = sns.heatmap(
            self.df.T,
            cmap=cmap,
            yticklabels=True,
            cbar_kws={'label': value_label},
            center=0,
            linewidths=0.01 if gaps else 0,
            linecolor='whitesmoke',
            # square=True
        )

        # Set all ticks
        ax.tick_params(
            left=True,
            bottom=True,
            right=True,
            top=True,
            labelright=True,
            labeltop=True,
            labelrotation=0
        )

        # Plot boxes for wild types
        if wts:
            wt_labels = dict(zip(
                zip(self.df.index, self.wts),
                self.n_residues * ['magenta']
            ))
            wt_labels = {k: v for k, v in wt_labels.items() if k[1] != '-'}
            self._plot_labels(ax, wt_labels)

        # Plot specified labels
        if labels is not None:
            self._plot_labels(ax, labels)

        # Set title
        ax.set_title(self.name)

        return ax

    def _plot_labels(self, ax, labels: Dict[tuple, str]):
        for (pos, letter), color in labels.items():
            letter_id = list(self.df.columns).index(letter)
            ax.add_patch(Rectangle(
                (pos, letter_id),
                1, 1,
                fill=False,
                edgecolor=color,
                lw=0.25)
            )
        return ax

    def plot_label_swarm(
        self,
        labels: Dict[tuple, str],
        ax: Any = None
    ):
        values = [self.df.loc[pos, mut] for pos, mut in labels.keys()]
        ax = sns.swarmplot(
            x=values,
            hue=labels.values(),
            palette=['green', 'red'],
            ax=ax
        )
        ax.axvline(x=0, color='black', linestyle='dashed')
        ax.legend(loc='center left', title='Label', bbox_to_anchor=(1, 0.5))
        ax.set_title(self.name)
        return ax

    def _test_df_format(self):
        assert self.df.columns.tolist() == AMINO_ACID_CODES_1

    def _ipython_display_(self):
        display(self.df)

    def __eq__(self, other):
        if type(other) is type(self):
            return self.df.equals(other.df) and self.name == other.name
        return False
    
    def operate(self, other, op, repr):
        #common logic here
        retval = copy.deepcopy(self)
        retval.df = op(self.df, other.df)
        retval.name = f'{self.name} {repr} {other.name}'
        return retval

    def __add__(self, other):
        return self.operate(other, operator.add, '+')

    def __sub__(self, other):       
        return self.operate(other, operator.sub, '-') 

    def copy(self):
        return copy.deepcopy(self)
    
    def to_csv(self, path: Union[Path, str]):
        if isinstance(path, str):
            path = Path(path)
        self.df.to_csv(path)


class NormalizedImprovementScoringMatrix(ScoringMatrix):
    SCORE_TITLE = 'Normalized Improvement Score'

    def __init__(
        self,
        df: pd.DataFrame,
        wts: str = None,
        name: str = '',
        weight: float = 1.0,
        normalized: bool = False,
        improvements: bool = False,
        *args,
        **kwargs
    ):
        super().__init__(df, wts, name, weight, *args, **kwargs)
        if not normalized:
            self.standardize(shift=False)
        if not improvements:
            self.subtract_wts()

    def plot(
        self,
        cmap: Any = sns.blend_palette(
            ['#df0101', '#f5f6ce', '#31b404'], as_cmap=True
        ),
        value_label: str = 'Normalized improvement over wild type',
        **kwargs
    ):
        return super().plot(
            cmap=cmap,
            value_label=value_label,
            **kwargs
        )


class RankingMatrix:
    def __init__(self):
        raise NotImplementedError()


class ScoringMatrixCollection:
    def __init__(
        self,
        matrices: List[ScoringMatrix]
    ):
        self.matrices = matrices
        self.matrix_dict = {
            matrix.name: matrix for matrix in matrices
        }
        self.group_dict = defaultdict(list)

        self._test_integrity()

    def groupby(self, matrix_name_to_group_name: Callable[[str], str]):
        # Group matrices according to provided mapping
        for matrix in self.matrices:
            group_name = matrix_name_to_group_name(matrix.name)
            self.group_dict[group_name].append(matrix.name)

    def average(
        self,
        weights: Iterable[float] = None,
        dropna: Literal['none', 'full_pos_or'] = 'full_pos_or'
    ):
        """
        :param weights:
        :param dropna: 'full_pos_or' means that the position is dropped if at
            least of one the matrices has this position full of NaNs.
        :return:
        """
        # Init weights
        if weights is None:
            weights = [matrix.weight for matrix in self.matrices]

        # Calculate average ignoring NaNs
        tensor = np.array(list(map(lambda m: m.df.to_numpy(), self.matrices)))
        masked_tensor = np.ma.masked_array(tensor, np.isnan(tensor))
        matrix = np.ma.average(masked_tensor, axis=0, weights=weights)
        matrix = matrix.filled(np.nan)

        # Restore format
        template = self.matrices[0]

        kwargs = dict()
        if isinstance(template, NormalizedImprovementScoringMatrix):
            kwargs |= dict(normalized=True, improvements=True)

        matrix = type(template)(
            pd.DataFrame(
                matrix,
                columns=template.df.columns,
                index=template.df.index
            ),
            wts=template.wts,
            weight=float(np.mean(weights)),
            **kwargs
        )

        # Handle NaNs
        if dropna == 'full_pos_or':
            pos_masks = [
                matrix.df.isna().sum(axis=1) == matrix.n_letters
                for matrix in self.matrices
            ]
            pos_mask = np.logical_or.reduce(pos_masks)
        elif dropna == 'none':
            pos_mask = None
        else:
            raise NotImplementedError()
        if pos_mask is not None:
            matrix.df.loc[pos_mask] = np.nan

        return matrix

    def average_groups(self, *args, **kwargs):
        averaged_groups = []
        for group_name in self.group_names():
            averaged = self.group(group_name).average(*args, **kwargs)
            averaged.name = group_name
            averaged_groups.append(averaged)
        averaged_groups = ScoringMatrixCollection(averaged_groups)
        return averaged_groups

    def fit_weights(
        self,
        labels: Dict[tuple, Union[int, float]],
        balanced: bool = False,
        n_epochs: int = 100
    ):
        mutations = labels.keys()
        y = list(labels.values())
        y = torch.tensor(y, dtype=torch.double)

        # Construct loss weights
        loss_weights = y.mean()*(1-y) + (1-y.mean())*y

        # Convert weights to Torch parameters
        weights = [matrix.weight for matrix in self.matrices]
        weights = torch.tensor(weights)
        weights.requires_grad_()

        # Optimize weights
        optimizer = torch.optim.Adam([weights])
        loss = torch.nn.BCELoss(weight=loss_weights if balanced else None)
        for i in range(n_epochs):
            optimizer.zero_grad()
            weights_normalized = torch.softmax(weights, dim=-1)

            # Predict
            y_pred = []
            for matrix, weight in zip(self.matrices, weights_normalized):
                matrix_pred = matrix.predict(mutations)
                matrix_pred = torch.tensor(matrix_pred)
                matrix_pred *= weight
                y_pred.append(matrix_pred)
            y_pred = torch.stack(y_pred).mean(0)
            # print(torch.mean((torch.tensor(y_pred >= 0).float() == y).float()))
            y_pred = torch.sigmoid(y_pred)

            # Optimize
            output = loss(y_pred, y)
            output.backward()
            optimizer.step()

        # Normalize weights
        weights = torch.softmax(weights, dim=-1)

        # Set updated weights
        for matrix, weight in zip(self.matrices, weights):
            matrix.weight = weight.detach().numpy().item()

        return self

    def correlate(
        self,
        positions: Collection[int] = None,
        plot: bool = True
    ):
        # Set to all positions if not specified
        if positions is None:
            positions = self.matrices[0].df.index.tolist()

        # Calculate pair-wise correlation dataframe over positions
        df_corrs = []
        for m0 in self.matrices:
            for m1 in self.matrices:
                df0, df1 = m0.df, m1.df
                df0 = df0.loc[positions]
                df1 = df1.loc[positions]
                corrs_pair = []
                for (i, r0), (_, r1) in zip(df0.iterrows(), df1.iterrows()):
                    res = stats.spearmanr(r0, r1, nan_policy='omit')
                    corrs_pair.append(res.correlation)
                df_corrs.append([m0.name, m1.name] + corrs_pair)
        df_corrs = pd.DataFrame(
            df_corrs,
            columns=['Matrix 1', 'Matrix 2'] + positions
        )

        # Calculate mean pair-wise correlation over positions
        df_corrs['Mean Spearman Correlation'] =\
            df_corrs.apply(lambda x: x[positions].mean(), axis=1)

        # Plot
        if plot:
            df_corrs_melt = pd.melt(df_corrs, id_vars=['Matrix 1', 'Matrix 2'])
            df_corrs_crosstab = pd.crosstab(
                df_corrs['Matrix 1'],
                df_corrs['Matrix 2'],
                df_corrs['Mean Spearman Correlation'],
                aggfunc='mean'
            )
            # TODO Enhance plotting
            # TODO Remove
            idx = sorted(list(map(lambda x: x[::-1] , df_corrs_crosstab.index)))
            idx = list(map(lambda x: x[::-1] , idx))
            cmap = sns.color_palette("blend:#df0101,#f5f6ce,#31b404", as_cmap=True)
            df_corrs_crosstab.loc[idx, idx]
            fig = sns.heatmap(
                df_corrs_crosstab.loc[idx, idx],
                annot=True,
                square=True,
                cbar_kws={"shrink": 0.75},
                annot_kws={"size": 12},
                cmap=cmap,
                center=0
            )
            plt.gca().tick_params(labelsize=12)
            cbar = plt.gca().collections[0].colorbar
            cbar.ax.tick_params(labelsize=12)
            plt.gca().tick_params(
                left=False,
                bottom=False,
                right=False,
                top=False,
                labelright=False,
                labeltop=True,
                labelleft=True,
                labelbottom=False,
                # labelrotation=90
            )
            plt.xticks(rotation=90) 
            plt.xlabel(None)
            plt.ylabel(None)

            return df_corrs, fig

        return df_corrs


    def names(self):
        return [matrix.name for matrix in self.matrices]

    def group_names(self):
        return list(self.group_dict.keys())

    def group(self, group_name: str):
        return self[self.group_dict[group_name]]

    def plot(self, *args, **kwargs):
        for matrix in self.matrices:
            matrix.plot(*args, **kwargs)

    def _test_integrity(self):
        # Test collection is not empty
        assert len(self.matrices)

        # Test column names in each matrix are identical
        columns = self.matrices[0].df.columns
        for matrix in self.matrices[1:]:
            assert all(columns == matrix.df.columns)

    def __getitem__(self, item):
        if isinstance(item, int):
            return self.matrices[item]
        if isinstance(item, list):
            return ScoringMatrixCollection([self[i] for i in item])
        if isinstance(item, str):
            return self.matrix_dict[item]
        raise ValueError('Wrong `item` value.')

    # TODO Group name?
    def __setitem__(self, key, value):
        if isinstance(key, int):
            self.matrices[key] = value
            self.matrix_dict[self.names()[key]] = value
        elif isinstance(key, str):
            self.matrix_dict[key] = value
            if key in self.names():
                self.matrices[self.names().index(key)] = value
            else:
                self.matrices.append(value)
        else:
            raise ValueError('Wrong `item` value.')

    def append(self, matrix, group_name=None):
        self.matrices.append(matrix)
        self.matrix_dict[matrix.name] = matrix
        if group_name is not None:
            self.group_dict[group_name].append(matrix.name)
        self._test_integrity()

    def __delitem__(self, key):
        self.matrices = type(self.matrices)(
            [m for m in self.matrices if m.name != key]
        )
        del self.matrix_dict[key]
        for group_name in self.group_dict:
            if key in self.group_dict[group_name]:
                self.group_dict[group_name].remove(key)

    def __len__(self):
        return len(self.matrices)
