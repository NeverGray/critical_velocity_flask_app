# Project Name: Beyer-Stacey=Brenn CV Calculator
# Description: A web application to calculate the critical velocity based on Beryer-Stacey-Brenn
# Copyright (c) 2025 Justin Edenbaum, Never Gray
#
# This file is licensed under the MIT License.
# You may obtain a copy of the license at https://opensource.org/licenses/MIT

from flask import Flask, render_template, request, redirect, url_for, flash
import io
import base64
import matplotlib
matplotlib.use('Agg')  # Use non-GUI backend

from src.critical_velocity import Fire, Tunnel, iterate_critical_velocity, fire_response

app = Flask(__name__)
app.secret_key = 'your_secret_key'  # Replace with a secure key in production

@app.route('/', methods=['GET', 'POST'])
def index():
    # Default values
    DEFAULTS = {
        'intensity': 2.25,      # Fire intensity (WW/m²)
        'width': 2.5,             # Max width of fire source (m)
        'epsilon': 0.2,          # Heat reduction due to imperfect combustion
        'eta': 0.0,               # Reduction of convective heat release
        'ambient_temp': 294,   # Temperature of upstream air (K)
        'ambient_density': 1.2    # Ambient density (kg/m³)
    }

    form_values = {}
    if request.method == 'POST':
        try:
            # User Input
            hrr = request.form['hrr']
            height = request.form['height']
            area = request.form['area']

            # Handle defaults for each variable
            use_default_intensity = 'use_default_intensity' in request.form
            intensity = DEFAULTS['intensity'] if use_default_intensity else request.form['intensity']

            use_default_width = 'use_default_width' in request.form
            width = DEFAULTS['width'] if use_default_width else request.form['width']

            use_default_epsilon = 'use_default_epsilon' in request.form
            epsilon = DEFAULTS['epsilon'] if use_default_epsilon else request.form['epsilon']

            use_default_eta = 'use_default_eta' in request.form
            eta = DEFAULTS['eta'] if use_default_eta else request.form['eta']

            use_default_ambient_temp = 'use_default_ambient_temp' in request.form
            ambient_temp = DEFAULTS['ambient_temp'] if use_default_ambient_temp else request.form['ambient_temp']

            use_default_ambient_density = 'use_default_ambient_density' in request.form
            ambient_density = DEFAULTS['ambient_density'] if use_default_ambient_density else request.form['ambient_density']

            hydraulic_diameter = request.form['hydraulic_diameter']

            # Save for repopulation
            form_values = {
                'height': height,
                'area': area,
                'hydraulic_diameter': hydraulic_diameter,
                'hrr': hrr,
                'intensity': intensity,
                'width': width,
                'epsilon': epsilon,
                'eta': eta,
                'ambient_temp': ambient_temp,
                'ambient_density': ambient_density,
                'use_default_intensity': use_default_intensity,
                'use_default_width': use_default_width,
                'use_default_epsilon': use_default_epsilon,
                'use_default_eta': use_default_eta,
                'use_default_ambient_temp': use_default_ambient_temp,
                'use_default_ambient_density': use_default_ambient_density
            }

            # Convert to float for calculation
            #Tunnel Properties
            height = float(height)
            area = float(area)
            hydraulic_diameter = float(hydraulic_diameter)
            #Fire Properties
            hrr = float(hrr)*1e6 # Convert MW to W
            intensity = float(intensity)*1e6 # Convert MW/m² to W/m²
            width = float(width)
            epsilon = float(epsilon)
            eta = float(eta)
            #Ambient Properties
            ambient_temp = float(ambient_temp)
            ambient_density = float(ambient_density)
            # Default Constants
            K_g = 1.0
            tol = 1e-6
            fire = Fire(hrr, intensity, width, epsilon, eta)
            tunnel = Tunnel("User Tunnel", height, area, hydraulic_diameter)
            critical_velocity = 2.0  # Initial guess
            critical_velocity, dt = iterate_critical_velocity(
                fire, tunnel, critical_velocity, ambient_temp, ambient_density, epsilon, eta, K_g, tol
            )

            # Format the results
            critical_velocity = round(critical_velocity, 3)
            dt = round(dt, 1)

            # Plotting
            plot_url = None
            try:
                fig = fire_response(
                    [tunnel],
                    1e6,  # min_hrr
                    200e6,  # max_hrr
                    1e-5,  # tol
                    float(ambient_temp),
                    101325.0,  # ref_pressure
                    float(intensity),
                    float(width),
                    float(epsilon),
                    float(eta),
                    1.0,  # K_g
                    for_web=True
                    )
                buf = io.BytesIO()
                fig.savefig(buf, format='png')
                buf.seek(0)
                plot_url = base64.b64encode(buf.getvalue()).decode('utf8')
            except Exception as e:
                plot_url = None

            return render_template(
                'index.html',
                critical_velocity=critical_velocity,
                delta_t=dt,
                form_values=form_values,
                DEFAULTS=DEFAULTS,
                plot_url=plot_url
            )

        except ValueError:
            flash('Please enter valid numerical values for all fields.')
            return redirect(url_for('index'))

    # For GET: set all use_default_* to True
    form_values = {
        'use_default_intensity': True,
        'use_default_width': True,
        'use_default_epsilon': True,
        'use_default_eta': True,
        'use_default_ambient_temp': True,
        'use_default_ambient_density': True,
        # Optionally, set the default values for the fields as well:
        'intensity': DEFAULTS['intensity'],
        'width': DEFAULTS['width'],
        'epsilon': DEFAULTS['epsilon'],
        'eta': DEFAULTS['eta'],
        'ambient_temp': DEFAULTS['ambient_temp'],
        'ambient_density': DEFAULTS['ambient_density'],
    }
    return render_template('index.html', form_values=form_values, DEFAULTS=DEFAULTS)

if __name__ == '__main__':
    app.run(debug=False)