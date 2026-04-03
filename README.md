# Sheet Pile Design Streamlit App

A Streamlit prototype for steel sheet pile / soldier pile retaining wall design based on the logic shown in the provided screenshots.

## Included features
- Input panel similar to the sample sheet
- Force calculation at cantilever bottom
- Section limitation checks
- Approximate AISC 360 H1 ASD interaction check
- Iterative embedment depth calculation following the shown equation format
- Built-in mini section database plus manual section property mode
- CSV export of the current design result

## Files
- `app.py` - main Streamlit application
- `requirements.txt` - Python dependencies
- `data/section_database.csv` - editable sample steel section library

## Run locally
```bash
pip install -r requirements.txt
streamlit run app.py
```

## Notes
- This version follows the calculation flow visible in the screenshots.
- The restrained top force inputs are included in the UI, but are not yet directly used in the core force equations because that branch was not fully visible in the screenshots.
- The AISC strength check is implemented as a practical first-pass design aid. Final design should still be verified with full code checks and project-specific geotechnical assumptions.
