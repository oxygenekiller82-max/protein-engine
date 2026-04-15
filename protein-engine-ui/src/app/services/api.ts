import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { Observable } from 'rxjs';

@Injectable({ providedIn: 'root' })
export class ApiService {
  private apiUrl = 'http://localhost:5000/api/start_search'; 

  constructor(private http: HttpClient) {}

  startSearch(data: any): Observable<any> {
    return this.http.post<any>(this.apiUrl, data);
  }
}