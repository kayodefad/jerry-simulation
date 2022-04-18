from flask import Blueprint, request, jsonify
from src.mold_func import get_temp_varies

api = Blueprint("api", __name__, url_prefix="/api/v1/")

@api.post('/img_path')
def simulation():
    heat_conductivity = request.json['heat_conductivity']
    initial_temperature = request.json['initial_temperature']
    length = request.json['length']
    breadth = request.json['breadth']
    height = request.json['height']
    number_of_hours = request.json['number_of_hours']

    imageUrl = get_temp_varies(
            Heat_conductivity_of_walls=float(heat_conductivity),
            Initial_Temperature=int(initial_temperature),
            Length=int(length),
            Breadth=int(breadth),
            Height=int(height),
            Number_of_hours_to_run_simulation=int(number_of_hours),
        )

    return jsonify({
        'imageUrl': imageUrl,
    }), 201

