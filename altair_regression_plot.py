

import altair as alt
import pandas as pd
import numpy as np
import statsmodels.api as sm

def generate_demo_data():
    # Generate some sample data
    np.random.seed(0)
    x_obs = np.linspace(0, 10, 100)
    y_obs = x_obs + np.random.normal(0, 1, 100)

    # Create a pandas DataFrame
    data = pd.DataFrame({'experimental': x_obs, 'predicted': y_obs})
    return data

def altair_regression_plot(data, x_col, y_col):
    # Add a constant for the regression model
    X = sm.add_constant(data[x_col])
    
    # Fit the OLS model
    model = sm.OLS(data[y_col], X).fit()

    # Get the predictions and confidence intervals
    predictions = model.get_prediction(X).summary_frame(alpha=0.05)
    
    # Add the predictions and CI bounds to the DataFrame
    data['fit'] = predictions['mean']
    data['ci_lower'] = predictions['mean_ci_lower']
    data['ci_upper'] = predictions['mean_ci_upper']

    # The scatter plot of the original data
    scatter = alt.Chart(data).mark_point().encode(
        x=x_col,
        y=y_col
    )

    # The regression line
    line = alt.Chart(data).mark_line(color='red').encode(
        x=x_col,
        y='fit'
    )

    # The confidence interval as a shaded area
    band = alt.Chart(data).mark_area(opacity=0.3).encode(
        x=x_col,
        y='ci_lower',
        y2='ci_upper'
    )

    # Layer the charts
    chart = scatter + band + line
    return chart
    
    # Save the chart to an HTML file
    #chart.save('regression_plot.html')
    
    #print("Regression plot saved to regression_plot.html")

if __name__ == "__main__":
    demo_df = generate_demo_data()
    altair_regression_plot(demo_df,"experimental","predicted")
    
    
