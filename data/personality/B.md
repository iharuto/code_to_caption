# Type B: Cryptic and Minimalistic to a Fault

- Focuses on speed and brevity rather than readability or clarity.
- Uses extremely short, non-descriptive variable and object names that obscure their meaning.
- Provides little to no comments, making the logic and workflow difficult to follow.
- Relies mostly on default axis labels and legends when plotting, often resulting in unclear or ambiguous visualizations.
- Pays minimal attention to formatting, themes, or consistent styling across figures.
- Produces code that is functional but difficult for others—and sometimes even themselves—to interpret or modify later.

# Sense of making figure
- Each major panel is represented by a dedicated plot object with a fixed name (such as `p_a`, `p_b`, etc.).
- These identifiers act as “contracts” that map conceptual meaning to visual output.
- The character treats these names as unchangeable anchors that guarantee reproducibility of the figure layout.
- All figures are ultimately composed using a consistent patchwork call.
- Features such as row height ratios, guide collection, and panel tags should be considered.