"""
ana.fanal_display
=================
Utilities to display notation-to-code mapping tables in Jupyter notebooks.

Helps students verify their computed values before writing them to collpars.py.
"""

from IPython.display import display, HTML


def display_collpars(entries, title="Variables to write to collpars.py"):
    """Render a notation-to-code table with the student's computed values.

    Parameters
    ----------
    entries : list of tuples
        Each element is ``(latex, varname, value)`` or
        ``(latex, varname, value, fmt)`` where

        - *latex*   – LaTeX math expression (without ``$``), e.g. ``r'\\epsilon^\\mathrm{Bi}_E'``
        - *varname* – Python variable name as it appears in collpars.py
        - *value*   – the computed value (number, string, tuple …)
        - *fmt*     – optional format string (default ``'.3g'`` for numbers)
    title : str, optional
        Table caption.
    """
    rows = []
    for entry in entries:
        if len(entry) == 4:
            latex, varname, value, fmt = entry
        else:
            latex, varname, value = entry
            fmt = None

        cell_value = _format_value(value, fmt)
        rows.append(f"<tr><td>${latex}$</td>"
                    f"<td><code>{varname}</code></td>"
                    f"<td style='text-align:right'>{cell_value}</td></tr>")

    html = (
        f"<table style='border-collapse:collapse; margin:12px 0'>"
        f"<caption style='font-weight:bold; text-align:left; "
        f"margin-bottom:6px'>{title}</caption>"
        f"<thead><tr>"
        f"<th style='border:1px solid #888; padding:4px 10px'>Math</th>"
        f"<th style='border:1px solid #888; padding:4px 10px'>Python variable</th>"
        f"<th style='border:1px solid #888; padding:4px 10px'>Your value</th>"
        f"</tr></thead><tbody>"
        + "\n".join(rows)
        + "</tbody></table>"
    )
    # apply common cell style
    html = html.replace("<td>", "<td style='border:1px solid #888; padding:4px 10px'>")
    display(HTML(html))


def _format_value(value, fmt=None):
    """Format a single value for display."""
    if value is None:
        return "&mdash;"
    if isinstance(value, str):
        return value
    if isinstance(value, tuple):
        if fmt:
            inner = ", ".join(f"{v:{fmt}}" for v in value)
        else:
            inner = ", ".join(f"{v}" for v in value)
        return f"({inner})"
    if isinstance(value, float):
        fmt = fmt or ".4g"
        return f"{value:{fmt}}"
    return str(value)
