from flask import Flask, render_template, request, redirect, url_for
from aspirin_main import aspirin_main

app = Flask(__name__)

# def hello_world():
#     return 'Hello World!'
@app.route('/')
def home():
    return render_template('index.html')

@app.route('/calculate', methods=['POST'])
def index():
    if request.method == 'POST':
        body_weight = float(request.form['body_weight'])
        dose_step = float(request.form['dose_step'])
        dose, time_next_dose, plot_data = aspirin_main(body_weight, dose_step)

        # Concatenate t1 and t2, y1 and y2 in Python
        plot_data['t'] = plot_data['t1'] + plot_data['t2']
        plot_data['y'] = plot_data['y1'] + plot_data['y2']

        return render_template('result.html', dose=dose, time_next_dose=time_next_dose, plot_data=plot_data)


if __name__ == '__main__':
    app.run()
