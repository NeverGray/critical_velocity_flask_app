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
import matplotlib.pyplot as plt

from src.critical_velocity import Fire, Tunnel, iterate_critical_velocity, plot_critical_velocity

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
        'ambient_pressure': 101325    # Ambient pressure (Pa)
    }
    
    # Form field definitions - this makes the HTML much more maintainable
    FORM_FIELDS = [
        {
            'category': 'Tunnel Properties',
            'fields': [
                {
                    'name': 'height',
                    'symbol': r'\( L_n \)',
                    'description': 'Characteristic length (height) for natural convection (for pool fires: tunnel height from the base of the fire to highest point at the ceiling)',
                    'unit': 'm',
                    'has_default': False,
                    'required': True,
                    'greater_than_zero': True
                },
                {
                    'name': 'area',
                    'symbol': r'\( A \)',
                    'description': 'Tunnel cross-sectional area at the fire site',
                    'unit': 'm²',
                    'has_default': False,
                    'required': True,
                    'greater_than_zero': True
                },
                {
                    'name': 'hydraulic_diameter',
                    'symbol': r'\( D_h \)',
                    'description': 'Tunnel hydraulic diameter at the fire location',
                    'unit': 'm',
                    'has_default': False,
                    'required': True,
                    'greater_than_zero': True
                }
            ]
        },
        {
            'category': 'Fire Properties',
            'fields': [
                {
                    'name': 'hrr',
                    'symbol': r'\( \dot{Q} \)',
                    'description': 'Total fire heat release rate',
                    'unit': 'MW',
                    'has_default': False,
                    'required': True,
                    'greater_than_zero': True
                },
                {
                    'name': 'intensity',
                    'symbol': r'\( I_{fire} \)',
                    'description': 'Fire intensity (heat release rate per unit area)<br><small>e.g. pool fire with oil fuel 2.25 MW/m²</small>',
                    'unit': 'MW/m²',
                    'has_default': True,
                    'default_value': DEFAULTS['intensity'],
                    'required': True
                },
                {
                    'name': 'width',
                    'symbol': r'\( W_{fire} \)',
                    'description': 'Maximum width of the fire source<br><small>e.g. fire pans of Memorial Tunnel tests were 3.66 m wide</small>',
                    'unit': 'm',
                    'has_default': True,
                    'default_value': DEFAULTS['width'],
                    'required': True
                },
                {
                    'name': 'epsilon',
                    'symbol': r'\( \epsilon \)',
                    'description': 'Heat reduction due to imperfect combustion and thermal radiation',
                    'unit': '-',
                    'has_default': True,
                    'default_value': DEFAULTS['epsilon'],
                    'required': True
                },
                {
                    'name': 'eta',
                    'symbol': r'\( \eta \)',
                    'description': 'Reduction of the convective heat release rate due to water mist systems or sprinklers',
                    'unit': '-',
                    'has_default': True,
                    'default_value': DEFAULTS['eta'],
                    'required': True
                }
            ]
        },
        {
            'category': 'Ambient Properties',
            'fields': [
                {
                    'name': 'ambient_temp',
                    'symbol': r'\( T_a \)',
                    'description': 'Temperature of upstream air',
                    'unit': 'K',
                    'has_default': True,
                    'default_value': DEFAULTS['ambient_temp'],
                    'required': True,
                    'greater_than_zero': True
                },
                {
                    'name': 'ambient_pressure',
                    'symbol': r'\( p_a \)',
                    'description': 'Ambient Pressure',
                    'unit': 'Pa',
                    'has_default': True,
                    'default_value': DEFAULTS['ambient_pressure'],
                    'required': True
                }
            ]
        }
    ]

    form_values = {}
    if request.method == 'POST':
        try:
            # Process form data dynamically using our field definitions
            form_values = {}
            processed_values = {}
            
            # Get all field names from our structure
            all_fields = []
            for category in FORM_FIELDS:
                all_fields.extend(category['fields'])
            
            # Process each field
            for field in all_fields:
                field_name = field['name']
                use_default_key = f'use_default_{field_name}'
                
                # Check if using default
                use_default = use_default_key in request.form
                form_values[use_default_key] = use_default
                
                # Get the value
                if use_default and field.get('has_default'):
                    value = field['default_value']
                else:
                    value = request.form[field_name]
                
                form_values[field_name] = value
                processed_values[field_name] = value

            # Convert to float for calculation
            #Tunnel Properties
            height = float(processed_values['height'])
            area = float(processed_values['area'])
            hydraulic_diameter = float(processed_values['hydraulic_diameter'])
            #Fire Properties
            hrr = float(processed_values['hrr'])*1e6 # Convert MW to W
            intensity = float(processed_values['intensity'])*1e6 # Convert MW/m² to W/m²
            width = float(processed_values['width'])
            epsilon = float(processed_values['epsilon'])
            eta = float(processed_values['eta'])
            #Ambient Properties
            ambient_temp = float(processed_values['ambient_temp'])
            ambient_pressure = float(processed_values['ambient_pressure'])
            
            # Input Validation: Strictly positive values check
            if height <= 0 or area <= 0 or hydraulic_diameter <= 0 or ambient_temp <= 0:
                 flash('Height, Area, Hydraulic Diameter, and Ambient Temperature must be strictly positive.')
                 return render_template('index.html', form_values=form_values, DEFAULTS=DEFAULTS, FORM_FIELDS=FORM_FIELDS)

            fire = Fire(hrr, intensity, width, epsilon, eta)
            tunnel = Tunnel("Web Tunnel", height, area, hydraulic_diameter)

            # Plotting
            plot_url = None
            critical_velocity = 0
            fire_hrr = 0
            oxygen_depletion_msg = ""
            converging_msg = ""

            try:
                fig, fire_hrr, critical_velocity, oxygen_depletion_msg, converging_msg = plot_critical_velocity(
                    tunnel,
                    fire,
                    ambient_temp=ambient_temp,
                    ambient_pressure=ambient_pressure,
                )
                buf = io.BytesIO()
                fig.savefig(buf, format='png')
                buf.seek(0)
                plot_url = base64.b64encode(buf.getvalue()).decode('utf8')
                plt.close(fig)

            except Exception as e:
                plot_url = None
                flash(f"An error occurred during calculation: {str(e)}")
                return render_template('index.html', form_values=form_values, DEFAULTS=DEFAULTS, FORM_FIELDS=FORM_FIELDS)

            return render_template(
                'index.html',
                critical_velocity=round(critical_velocity, 3),
                fire_hrr=round(fire_hrr / 1e6, 1),  # Convert back to MW for display
                form_values=form_values,
                DEFAULTS=DEFAULTS,
                FORM_FIELDS=FORM_FIELDS,
                plot_url=plot_url,
                oxygen_depletion_msg=oxygen_depletion_msg
            )

        except ValueError:
            flash('Please enter valid numerical values for all fields.')
            return render_template('index.html', form_values=form_values, DEFAULTS=DEFAULTS, FORM_FIELDS=FORM_FIELDS)

    # For GET: set defaults dynamically based on field definitions
    form_values = {}
    for category in FORM_FIELDS:
        for field in category['fields']:
            field_name = field['name']
            if field.get('has_default'):
                form_values[f'use_default_{field_name}'] = True
                form_values[field_name] = field['default_value']
    
    return render_template('index.html', form_values=form_values, DEFAULTS=DEFAULTS, FORM_FIELDS=FORM_FIELDS)

if __name__ == '__main__':
    app.run(debug=False)