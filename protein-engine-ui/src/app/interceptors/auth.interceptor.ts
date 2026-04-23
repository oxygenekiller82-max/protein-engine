import { HttpInterceptorFn,HttpErrorResponse } from '@angular/common/http';

import { inject } from '@angular/core';

import { catchError, switchMap, throwError } from 'rxjs';

import { Auth } from '../services/auth';

export const authInterceptor: HttpInterceptorFn = (req, next) => {
    const authService = inject(Auth);
    const token = localStorage.getItem('accessToken');

    let request = req;

    if (token) {
      const payload = JSON.parse(atob(token.split('.')[1]));
      const expiresAt = new Date(payload.exp * 1000);
      const now = new Date();
      console.log(`[JWT] request: ${req.url}`);
      console.log(`[JWT] token expires: ${expiresAt}`);
      console.log(`[JWT] now: ${now}`);
      console.log(`[JWT] expired: ${now > expiresAt}`);
    } else {
      console.warn(`[JWT] NO TOKEN for request: ${req.url}`);
    }

    //have token ?-> clone the request -> auth bearer token 
    if (token) {
        request = req.clone({
          setHeaders: {
            Authorization: `Bearer ${token}`
          },
          body: req.body
        });
    }
    
    if (req.url.includes('/refresh')) {
      return next(request);
    }

    //401 -> refresh that token sir 
    return next(request).pipe(
      catchError((error) => {
          if (error instanceof HttpErrorResponse && error.status === 401) {
              const refreshToken = localStorage.getItem('refreshToken');

              if (refreshToken) {
                //get new access token
                return authService.refresh(refreshToken).pipe(
                    switchMap((res: any) => {
                      localStorage.setItem('accessToken', res.accessToken);
                      //RETRY ORIGNAL REQUEST but with new token!
                      const retryReq = req.clone({
                        setHeaders: { Authorization: `Bearer ${res.accessToken}` }
                    });

                    return next(retryReq);
                  }),
                  catchError((err) => {
                    // Refresh token failed/expired too? Logout.
                    authService.logout();
                    return throwError(() => err);
                })
              );
            }
          }
          return throwError(() => error);
        })
      );
    };



