# Critical Velocity Flask Application

This project is a Flask web application that calculates the critical velocity of a fire in a tunnel based on user input. The calculatioin is in the `src/critical_velocity.py` file.

## Project Structure

```
critical_velocity_flask_app
├── app.py                  # Main entry point of the Flask application (Directions for the Web-site)
├── requirements.txt        # Lists the dependencies required for the project
├── src
│   ├── critical_velocity.py # Contains the logic for critical velocity calculation
├── templates
│   └── index.html          # HTML template for the web page
└── README.md               # Documentation for the project
```

## Installation

1. Clone the repository:
   ```
   git clone <repository-url>
   cd critical_velocity_flask_app
   ```

2. Create a virtual environment (optional but recommended):
   ```
   python -m venv venv
   source venv/bin/activate  # On Windows use `venv\Scripts\activate`
   ```

3. Install the required dependencies:
   ```
   pip install -r requirements.txt
   ```

## Usage

1. Run the Flask application:
   ```
   python app.py
   ```

2. Open your web browser and go to `http://127.0.0.1:5000`.

3. Fill in the form with the required parameters such as fire heat release rate, tunnel height, area, and other necessary values.

4. Submit the form to receive the calculated critical velocity.

## Contributing

Contributions are welcome! Please feel free to submit a pull request or open an issue for any suggestions or improvements.

## License
Copyright (c) 2025 Justin Edenbaum, Never Gray
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.