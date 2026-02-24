import { Component, signal } from '@angular/core';
import { RouterOutlet } from '@angular/router';
import { ColorShifter } from "./color-shifter/color-shifter";

@Component({
  selector: 'app-root',
  imports: [RouterOutlet, ColorShifter],
  templateUrl: './app.html',
  styleUrl: './app.css'
})
export class App {
  protected readonly title = signal('AAPFA');
}
