import { HttpClient } from '@angular/common/http';
import { Injectable,inject } from '@angular/core';
import { Observable } from 'rxjs';
import { Router } from '@angular/router';
import { SequenceService } from '../services/sequence.service';

@Injectable({
  providedIn: 'root',
})

//auth service actually..
export class Auth {
  private http = inject(HttpClient);
  private router = inject(Router);
  private seqService = inject(SequenceService);

  //spring URl is what..
  private readonly API_URL = 'http://localhost:8081/api/auth';


  //now how does this work.. 
  //Angular -> Asynchornous web requests 
  //Observable -> when u cann the method, server hasn't reponded yet 
  //it will promise to give it when it arrives or fails 
  //how to watch out for it ? ->.subscribe in the component! 

  register(signupData: any): Observable<any> {
    return this.http.post(`${this.API_URL}/register`, signupData);
    //signUpData = JSON object -> name email password..
  }

  login(loginData: any): Observable<any> {
    return this.http.post(`${this.API_URL}/login`, loginData);
  }

  //login = clear it all + THE BEHAVIOR SUBJCET!!  -> the .clearData in the sequence service
  logout():void{
    localStorage.removeItem('accessToken');
    localStorage.removeItem('refreshToken');
    this.seqService.clearData();

    this.router.navigate(['/login']);

  }

  getToken(): string | null {
    return localStorage.getItem('token');
  }
  //helper, token exists ?
  isLoggedIn(): boolean {
    return !!this.getToken(); // returns true if token exists, false if null
  }

  refresh(refreshToken: string): Observable<any> {// spring /refresh end point yh
    return this.http.post(`${this.API_URL}/refresh`, { refreshToken });
  }



}
