import { Component, OnInit } from '@angular/core';
import { Router } from '@angular/router';

@Component({
  selector: 'app-dashboard',
  templateUrl: './dashboard.html',
  styleUrls: ['./dashboard.css']
})
export class DashboardComponent implements OnInit {
  userName: string = '';

  stats = {
    sequences: 24,
    avgTime: 2.4,
    successRate: 87,
    explored: 15420
  };

  recentSequences = [
    { code: 'Ala-Gly-Ser-Val-Lys', hydro: '-2.1', mass: '532.6', binding: '124.3', date: '2024-01-15' },
    { code: 'Leu-Pro-Asn-Glu-Arg', hydro: '-1.8', mass: '598.2', binding: '156.7', date: '2024-01-14' },
    { code: 'Phe-Tyr-Trp-Cys-Met', hydro: '1.2', mass: '712.4', binding: '98.2', date: '2024-01-13' }
  ];

  constructor(private router: Router) {}

  ngOnInit() {
    this.userName = localStorage.getItem('user') || 'Chercheur';
    if (!localStorage.getItem('token')) {
      this.router.navigate(['/auth']);
    }
  }

  logout() {
    localStorage.removeItem('token');
    localStorage.removeItem('user');
    this.router.navigate(['/auth']);
  }
}