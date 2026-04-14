#API Comms file! 
import sys 
import json 
from AA_Math import start_search

from flask import Flask, request, jsonify, Response

from AA_Validation import compare_sequences
from AA_Validation import getSequenceData
from AA_Validation import reverse_mapping

app=Flask(__name__)

#IF sueccess -> STERAMILE the potentially 499K dicts file via the API

def streamline_DFS_history(DFS_history,stats,message,sequence): 
    """IF sucuessful DFS can be JUST below the LIMIT for calls -> gigantic list!
    -> must steamline it and not send in one Gb + file"""
    yield( 
        '{\n' + #JSON PARSER IGNORE \n 
        ' "status": "done" ,\n' +
        ' "message": ' + json.dumps(message) + ',\n' +
        ' "sequence": ' + json.dumps(sequence) + ',\n' +
        ' "history": \n[' #new line after history maybe ?
        )
    
     #creates JSOn object! wo


    for i, frame in enumerate(DFS_history):
        #1 dict -> one stream 
        frame='    '+json.dumps(frame) #4spaces
        # comma ->all elements except first
        if i>0:
            yield ",\n" + frame 
        else:
            yield frame

    #close 
    yield '],"stats": ' +json.dumps(stats) + '\n' + '}'


@app.route('/generate',methods=['POST'])
def generate(): 
    try:
        #1-get input from Spring Boot eeh..
        #input_json=sys.argv
        #params=json.loads(input_json)
        #-> request now 
        data=request.json 

        length=data.get('target_length',100)
        biological_switch=data.get('biological_switch',True)

        user_targets=data.get('targets',data)

        targets = data.get('targets', {})

        #STREAMLINE if done
        result=start_search(targets, length, biological_switch)
        status= result.get("status")
        sequence=result.get("final_sequence")
        DFS_history=result.get("history")
        message=result.get("message")
        stats= result.get("stats")

        if status=="done":
            return Response(
                streamline_DFS_history(DFS_history,stats,message,sequence),
                mimetype='application/json', #= label for the data huh..
                #application/json=> even tho just sending small batches of data 
                #final RESULT is a json object else SPring will think it's a text or vid or
                headers={"Content-Disposition": 'attachment; filename="PFA_DFS_search_log.log"'}
                  #-> Postman downloads that file -> + fixes .log extension!
            )
        elif status in ["length_error","Configuration Error"]: 
            #LENGTH_EXCEEDED, "nothing_found"
            return Response(
                json.dumps({
                "status": status,
                "sequence": [],
                "history":DFS_history,
                "message":message,
                "stats":stats
            }, indent=4), #indent 4 -> new line maybe ? o.o
            mimetype='application/json',
            headers={"Content-Disposition": 'attachment; filename="PFA_DFS_search_log.log"'}
        )

        #3-print JSON to Java 
        #print(json.dumps(results))    
        else:
            return Response(
            json.dumps({
            "status": status,
            "sequence": [],
            "history": DFS_history,
            "message": message,
            "stats": stats
            }, indent=4) ,
            mimetype='application/json',
            headers={"Content-Disposition": 'attachment; filename="PFA_DFS_search_log.log"'}
            ) 
    
    except Exception as e: #Error handling
        #print(json.dumps({"error":str(e)}))
        print(f"Error occurred: {e}") 
        return jsonify({"error": str(e)}), 500


@app.route('/generate_arrays',methods=['POST'])
def get_comparsion_data(): 
    #JSON -> body data
    data=request.get_json()

    #3letter ARRAYS now
    target_seq=data.get('target_seq')
    generated_seq=data.get('generated_seq')

    #jusst checking...
    if not target_seq or not generated_seq:
        return jsonify({"error": "Missing sequences."}), 400
    
    try: 
        #target 3->1 letter 
        target_seq_final = "".join([reverse_mapping[AA] for AA in target_seq])
        #generated 3->1 letter 
        generated_seq_final = "".join([reverse_mapping[AA] for AA in generated_seq])

        results = compare_sequences(target_seq, generated_seq)

        return jsonify(results)
    
    except KeyError as e: 
        return jsonify({"error": f"Invalid Amino Acid: {str(e)}"}), 400
    
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route('/get_peptide_stats', methods=['POST'])
def get_peptide_stats():
    data = request.get_json('sequence') #post body 

    sequence=data.get('sequence')#ARRAY ! of 3 letter codes

    if not sequence:
        return jsonify({"error": "No sequence found"}), 400
    
    try:
        #reverse mapping 3 letter -> 1 letter 
        sequence_string="".join([reverse_mapping[AA] for AA in sequence ])

        stats = getSequenceData(sequence_string.upper())
        return jsonify(stats)
    
    except KeyError:
        return jsonify({"error": "Invalid Amino Acid code found"}), 400
    
    except Exception as e:
        return jsonify({"error": str(e)}), 500

if __name__=="__main__":
    #main()
    app.run(port=5000,debug=False,use_reloader=False)

#why the if main ? 
#-> import a file in python => IT RUNS all the code in it 
#with this __main__ -> won't just start running the DFS file lol



