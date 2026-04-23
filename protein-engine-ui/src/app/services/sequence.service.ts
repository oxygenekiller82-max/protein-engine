import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { Observable, BehaviorSubject } from 'rxjs';

@Injectable({
    providedIn: 'root'
})
//what does this root mean ?
//-> this service is a SINGLETON ! 
//ONLY one instance of it in the entire app , apparently more efficient

//now we need a behaviorSubejct to store that log file..
// + Observable componene to watch it 
export class SequenceService{
    private apiUrl= 'http://localhost:8081/api/PFA/generate'

    private resultSource = new BehaviorSubject<any>(null);
    //that's it ?

    currentResult$ = this.resultSource.asObservable(); //OBSERVE IT o:
    //actually.. asObservable = allows it to READ only 

    //injection de dependace, httpClient to talk 
    constructor(private http: HttpClient) {}


    //now we generate 
    //generate(payload: any) -> now Home will call it , 
    //this.http.post -> use a POST request 
    //why <any> ? the log file is kinda complex ..so any 
    generate(payload:any): void{
        this.http.post<any>(this.apiUrl,payload).subscribe({
            //next -> uplaods the data onto the browser ram
            next: (data) =>{
                //WHEN it arrives -> push into the behaviour subject!
                this.resultSource.next(data);
            },
            error:(err)=>console.error("Log Pipeline failed",err)
        })
    }

    //logout = CLEAR THE BEHAVIOUR SUBJCET!
    clearData() {
        this.resultSource.next(null);
      }
}