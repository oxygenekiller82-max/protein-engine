import { Component } from '@angular/core';
import { CommonModule } from '@angular/common';

@Component({
  selector: 'app-logs',
  standalone: true,
  imports: [CommonModule],
  templateUrl: './logs.html',
  styleUrls: ['./logs.css']
})
export class LogsComponent {
  logs = [
    { id: 1, date: '2026-04-18', sequence: 'Sequence_Alpha', action: 'Generated', status: 'Success' },
    { id: 2, date: '2026-04-17', sequence: 'Sequence_Beta', action: 'Compared', status: 'Success' },
    { id: 3, date: '2026-04-16', sequence: 'Sequence_Gamma', action: 'Generated', status: 'Failed' },
    { id: 4, date: '2026-04-15', sequence: 'Sequence_Delta', action: 'Compared', status: 'Success' },
  ];
}