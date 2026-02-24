import { ComponentFixture, TestBed } from '@angular/core/testing';

import { ColorShifter } from './color-shifter';

describe('ColorShifter', () => {
  let component: ColorShifter;
  let fixture: ComponentFixture<ColorShifter>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      imports: [ColorShifter]
    })
    .compileComponents();

    fixture = TestBed.createComponent(ColorShifter);
    component = fixture.componentInstance;
    await fixture.whenStable();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
